from parameters import *

import math
import matplotlib.pyplot as plt
import ppigrf
from scipy.optimize import root_scalar
from scipy.spatial.transform import Rotation
from typing import Final


# Convert angular elements to radians and the semimajor axis from km to m.
KEPLER_INITIAL = {
    key: (math.radians(value) if key in {"INCL", "RAAN", "ARGP", "ANOM"}
          else value * 1000 if key == "SEMI"
          else value)
    for key, value in INITIAL_KEPLER_ELEMENTS.items()
}


def main():
    total_steps_per_cycle, steps_per_duty_cycle = validate_time_steps()
    SimulationController(total_steps_per_cycle, steps_per_duty_cycle)


def validate_time_steps():
        total_steps = DT_CONTROLLER / DT_SIMULATOR
        if total_steps.is_integer():
            controller_period = total_steps * DT_SIMULATOR
        else:
            total_steps = round(total_steps)
            controller_period = total_steps * DT_SIMULATOR
            print(f"Adjusted DT_CONTROLLER (controller period) to: {controller_period} s")

        # Ensure dissipation time at least 0.25 s
        if (1 - DUTY_CYCLE) * DT_CONTROLLER < 0.025:
            print("Duty cycle increased to meet minimum magnetorquer dissipation time requirements.")
            duty_steps = total_steps - math.ceil(0.025 / DT_SIMULATOR)
            return int(total_steps), int(duty_steps)
        duty_steps = DUTY_CYCLE * DT_CONTROLLER / DT_SIMULATOR
        if not duty_steps.is_integer():
            duty_steps = round(duty_steps)
            print(f"Adjusted DUTY_CYCLE to: {duty_steps / controller_period * DT_SIMULATOR}")
        return int(total_steps), int(duty_steps)


class SimulationController:
    start_time = START_DATE
    '''
    Controller is the class responsible for managing the simulation engine and 
    classes OrbitalModel, IGRFModel, SatelliteDynamics and BDotController
    '''
    
    def __init__(self, total_steps_per_cycle, steps_per_duty_cycle):
        self.start_time = self.real_time = START_DATE
        self.Ts = total_steps_per_cycle
        self.Tdis = steps_per_duty_cycle
        self.orbit = OrbitalModel()
        self.igrf = IGRFModel(self.orbit)
        self.sat_model = SatelliteDynamics(self.orbit)
        self.bdot = BDotController(total_steps_per_cycle, steps_per_duty_cycle, self.orbit)
        self.angular_velocity_times = [0.0]
        self.angular_velocity_history = [self.sat_model.w.copy()]

        # Start the simulation
        self.stepchange()

    def is_simulation_ended(self):
        end_time: Final = START_DATE + datetime.timedelta(hours=MAX_SIM_HOURS)
        return self.real_time >= end_time
        
    def stepchange(self):
        total_duration = MAX_SIM_HOURS * 3600
        progress_interval = total_duration / 20
        next_progress = progress_interval
        print(
            f"Simulation started for {MAX_SIM_HOURS} hours; "
            "the angular-velocity plot will open when it finishes."
            "\nSimulation progress: 0%"
        )
        
        while not self.is_simulation_ended():
            time_passed = self.real_time - self.start_time
            seconds_passed = time_passed.total_seconds()
            self.orbit.update_all_orbital_parameters(seconds_passed)
            
            curr_date = self.real_time.date()
            B_eci = self.igrf.return_B_ECI(curr_date)
            
            # Converts B vector from ECI to satellite body frame
            B_rot_matrix = self.sat_model.get_DCM_matrix()
            B_body_true = B_rot_matrix @ B_eci
            
            B_body_noisy = self.sat_model.add_B_noise(B_body_true)

            commanded_moment = self.bdot.stepchange(B_body_noisy)
            net_torque = self.sat_model.compute_net_torque(commanded_moment, B_body_true)
            self.sat_model.update_states(net_torque, B_rot_matrix, seconds_passed, DT_SIMULATOR)

            self.real_time += datetime.timedelta(seconds=DT_SIMULATOR)
            elapsed_seconds = (self.real_time - self.start_time).total_seconds()
            self.angular_velocity_times.append(elapsed_seconds)
            self.angular_velocity_history.append(self.sat_model.w.copy())
            if elapsed_seconds >= next_progress:
                progress = min(100, 100 * elapsed_seconds / total_duration)
                print(f"Simulation progress: {progress:.0f}%")
                next_progress += progress_interval

        self.display_angular_velocity()

    def display_angular_velocity(self):
        times_hours = np.asarray(self.angular_velocity_times) / 3600
        angular_velocity = np.asarray(self.angular_velocity_history)
        figure, axes = plt.subplots()
        for axis, component in enumerate(("x", "y", "z")):
            axes.plot(times_hours, angular_velocity[:, axis], label=f"ω{component}")
        axes.set_title("Angular velocity over time")
        axes.set_xlabel("Time (hours)")
        axes.set_ylabel("Angular velocity (rad/s)")
        axes.grid(True)
        axes.legend()
        figure.tight_layout()
        plt.savefig('simulation_results.png')
        plt.show()


class BDotController:
    max_moment: Final = MAX__MAGNETORQUER_MOMENT
    
    def __init__(self, total_steps_per_cycle, steps_per_duty_cycle, orbital_model):
        if (
            self.max_moment.shape != (3,)
            or not np.isfinite(self.max_moment).all()
            or np.any(self.max_moment <= 0)
        ):
            raise ValueError(
                "MAX__MAGNETORQUER_MOMENT must contain three positive per-axis limits."
            )

        self.mode = -1 # initialising mode
        self.mode_elapsed = 0 # number of steps elapsed during the current mode
        self.mode_length = total_steps_per_cycle * 4 # initialising mode is 4 cycles long (to get 5 initial measurements)
        self.mode_lengths = np.array([total_steps_per_cycle - steps_per_duty_cycle - 1, 
                                      1, 
                                      steps_per_duty_cycle]) # steps lengths for modes 0, 1, 2
        self.t_on: np.ndarray = np.full(3, np.nan) # takes only zero or positive values

        self.B_array = np.full((5, 3), np.nan) # unfiltered measurements of B
        self.moment = np.zeros(( 3 ))
        self.min_inertia_tensor = np.min(INERTIA_MATRIX.diagonal()) # matrix must be diagonalised
        self.orbit = orbital_model
        self.total_steps_per_cycle = total_steps_per_cycle

    def increment_mode(self):
        '''
        Alternate between different modes in each control cycle period.
        modes = {
            "off (initialising)": -1,
            "off (dissipating)": 0,
            "off (measuring)": 1, 
            "on (actuating)": 2 }
        '''
        mode_complete = (
            self.mode_elapsed >= self.mode_length if self.mode == -1
            else self.mode_elapsed + 1 >= self.mode_length
        )
        if not mode_complete:
            self.mode_elapsed += 1
            return

        next_mode = (self.mode + 1) % 3
        self.mode = next_mode
        self.mode_elapsed = 0
        self.mode_length = self.mode_lengths[next_mode]
    
    def return_filtered_points(self):
        '''
        Initially planned to use a Savitzky–Golay filter (polynomial fit) or a 
        moving average, but later abandoned due to time constraints.
        
        '''
        return self.B_array # does nothing

    def finite_differentiate(self, B_arr):
        '''
        Finite differentiates to find B_dot using the 5-point stencil method,
        taking only the rate of change of unit B vectors to enable a 
        projection-based controller.
        '''

        # Calculate last 5 B unit vectors
        norms = np.linalg.norm(B_arr, axis=-1, keepdims=True)
        norms[norms == 0] = 1.  # prevent division by 0
        B_unit_arr = B_arr / norms

        coefficients = np.array([3, -16, 36, -48, 25])
        b_dot = np.sum(coefficients[:, np.newaxis] * B_unit_arr, axis=0) / (12 * DT_CONTROLLER)
        return b_dot

    def update_desired_moment(self, B_reading):
        '''
        Output the final commanded moment, combining filtering (not implemented)
        and finite differentiation functions with a projection-based B-dot control
        law.
        '''

        # update self.B_array with raw value
        self.B_array = np.roll(self.B_array, -1, axis=0)
        self.B_array[-1] = B_reading # updates latest B value in array

        # return filtered data points
        filtered_B_arr = self.return_filtered_points()
        # calculate B dot
        Bdot = self.finite_differentiate(filtered_B_arr)
        
        # Calculate desired moment
        ## calculate the orbital angular velocity (time-constant)
        omega = (self.orbit.GP / self.orbit.SEMI ** 3) ** 0.5
        ## calculate the scalar gain
        k = 2 * omega * (1 + np.sin(self.orbit.INCL)) * self.min_inertia_tensor
        ## BDot control law
        moment = -k / np.linalg.norm(filtered_B_arr[-1]) * Bdot
        return moment
    
    def stepchange(self, B_reading):
        while self.mode >= 0 and self.mode_length == 0:
            self.increment_mode()

        if self.mode == -1: # initialising mode
            if self.mode_elapsed % self.total_steps_per_cycle == 0:
                self.moment = self.update_desired_moment(B_reading)
            self.increment_mode()
            return np.zeros(( 3 ))
        if self.mode == 0: # off (dissipating)
            self.increment_mode()
            return np.zeros(( 3 ))
        if self.mode == 1: # off (measuring)
            self.moment = self.update_desired_moment(B_reading)
            abs_moment_ratio = np.absolute(self.moment / self.max_moment)
            duty_time = self.mode_lengths[2]
            self.t_on = duty_time * np.clip(abs_moment_ratio, a_min=0, a_max=1)
            self.increment_mode()
            return np.zeros(( 3 ))
        if self.mode == 2: # on (actuating)
            actual_moment = np.full(3, np.nan)
            for i in range(len(actual_moment)):
                if self.mode_elapsed < self.t_on[i]:
                    actual_moment[i] = np.sign(self.moment[i]) * self.max_moment[i]
                else:
                    actual_moment[i] = 0.
            self.increment_mode()
            return actual_moment


class IGRFModel:
    omega_E = 7.2921151467 * 10**(-5) # mean angular velocity of the Earth (rad/s)
    def __init__(self, orbital_model):
        self.orbit = orbital_model

    def return_B_ECI(self, date: datetime.date):
        '''
        Series of transformations to convert Kepler elements in orbital model
        to position vector in geodetic coordinates (lon, lat, height) needed by
        the IGRF model. Then passes through ppigrf to get ENU (east, north, up)
        coordinates of B field vector, which is converted to ECI frame.
        '''

        # get ecef position vector
        eci = self.Kepler_to_ECI()
        ecef = self.ECI_to_ECEF(eci)

        lon, lat, height = self.ECEF_to_geodetic(ecef)

        # input to IGRF model
        igrf_date = datetime.datetime.combine(date, datetime.time.min)
        B_components = ppigrf.igrf(lon, lat, height, igrf_date)
        B_enu = np.array([component[0] for component in B_components])
        B_enu /= 10**9 # nanotesla to tesla (SI unit)

        lon_rad, lat_rad = math.radians(lon), math.radians(lat)
        sin_lon, cos_lon = math.sin(lon_rad), math.cos(lon_rad)
        sin_lat, cos_lat = math.sin(lat_rad), math.cos(lat_rad)
        enu_to_ecef = np.array([
            [-sin_lon, -sin_lat * cos_lon, cos_lat * cos_lon],
            [cos_lon, -sin_lat * sin_lon, cos_lat * sin_lon],
            [0.0, cos_lat, sin_lat],
        ])
        B_ecef = enu_to_ecef @ B_enu
        B_eci = self.ECI_to_ECEF(B_ecef, direction="reverse")
        return B_eci

    def ECEF_to_geodetic(self, coordinates):
        '''
        Convert an ECEF position in metres to geodetic coordinates.
        Written with AI assistance.
        '''

        x, y, z = coordinates
        semi_major_axis = 6378137.0
        eccentricity_squared = 6.69437999014e-3
        longitude = math.atan2(y, x)
        horizontal_radius = math.hypot(x, y)

        if horizontal_radius < 1e-12:
            latitude = math.copysign(math.pi / 2, z)
            height = abs(z) - semi_major_axis * math.sqrt(
                1 - eccentricity_squared
            )
        else:
            latitude = math.atan2(
                z, horizontal_radius * (1 - eccentricity_squared)
            )
            for _ in range(8):
                prime_vertical_radius = semi_major_axis / math.sqrt(
                    1 - eccentricity_squared * math.sin(latitude) ** 2
                )
                height = (
                    horizontal_radius / math.cos(latitude)
                    - prime_vertical_radius
                )
                latitude = math.atan2(
                    z,
                    horizontal_radius
                    * (
                        1
                        - eccentricity_squared
                        * prime_vertical_radius
                        / (prime_vertical_radius + height)
                    ),
                )

        return math.degrees(longitude), math.degrees(latitude), height / 1000
            
    def Kepler_to_ECI(self):
            '''
            Accesses instance attributes of Orbit and uses them to calculate 
            a position vector of the satellite in the perifocal frame, then 
            transform into the ECI frame (see logbook entry 2/10).
            '''

            r_peri = np.array([math.cos(self.orbit.ANOM), math.sin(self.orbit.ANOM), 0])
            r_peri *= self.orbit.SEMI * (1 - self.orbit.ECCT ** 2) / (1 + self.orbit.ECCT * math.cos(self.orbit.ANOM))
            rot = Rotation.from_euler("ZXZ", [self.orbit.RAAN, self.orbit.INCL, self.orbit.ARGP])
            r_ECI = rot.as_matrix() @ r_peri
            return r_ECI

    def ECI_to_ECEF(self, coordinates, direction="forward"):
        '''
        eta represents the change in rotation from the time ECI and ECEF frames coincided 
        (see entry 2/10 for details). direction="forward" represents ECI -> ECEF  
        transformation while direction="reverse" represents ECEF -> ECI.
        '''
        eta = self.omega_E * self.orbit.last_periapsis
        rot = Rotation.from_euler("Z", -eta)
        if direction == "forward":
            return rot.as_matrix() @ coordinates
        else:
            return rot.as_matrix().T @ coordinates


class OrbitalModel:
    J2_const = 1.08263 * 10**(-3)
    RE: Final = 6371 * 10**3            # radius of Earth
    GP: Final = 3.986004418 * 10**14    # gravitational parameter of Earth
    
    def __init__(self):
        self.SEMI = self.INCL = self.RAAN = self.ARGP = self.ANOM = self.ECCT = np.nan
        for key, value in KEPLER_INITIAL.items():
            setattr(self, key, value)

        self.eccentric_anomaly = 0.
        self.last_periapsis = 0.  # time since passed periapsis (coincides with ANOM = 0)
        self.dQ = np.nan # rate of change of RAAN (constant)
        self.dA = np.nan # rate of change of argument of periapsis (constant)
    
    def model_J2_effect(self):            
        if np.isnan(self.dQ):
            raan_numer = -3 * self.J2_const * self.RE**2 * self.GP**0.5 * math.cos(self.INCL)
            raan_denom = 2 * self.SEMI**3.5 * (1 - self.ECCT**2)**2
            self.dQ = raan_numer / raan_denom
        if np.isnan(self.dA):
            arg_perigee_numer = 3 * self.J2_const * self.RE**2 * self.GP**0.5 * (4 - 5 * math.sin(self.INCL)**2)
            arg_perigee_denom = 4 * self.SEMI**3.5 * (1 - self.ECCT**2)**2
            self.dA = arg_perigee_numer / arg_perigee_denom

        self.RAAN += DT_SIMULATOR * self.dQ
        self.RAAN = self.RAAN % (2 * math.pi)
        self.ARGP += DT_SIMULATOR * self.dA
        self.ARGP = self.ARGP % (2 * math.pi)

    def update_true_anomaly(self):
        mean_motion = (self.GP / self.SEMI**3) ** 0.5                           # units: rad/s
        mean_anomaly = (mean_motion * self.last_periapsis) % (2 * math.pi)      # units: rad

        def f(x):
            return (x - self.ECCT * math.sin(x) - mean_anomaly)
        
        self.eccentric_anomaly = root_scalar(f, x0=self.eccentric_anomaly, method="newton").root
        true_anomaly = 2 * math.atan2(
            math.sqrt(1 + self.ECCT)
            * math.sin(self.eccentric_anomaly / 2),
            math.sqrt(1 - self.ECCT)
            * math.cos(self.eccentric_anomaly / 2),
        )
        self.ANOM = true_anomaly % (2 * math.pi) # bounds true anomaly to [0, 2pi]

    def update_all_orbital_parameters(self, seconds_passed):
        self.last_periapsis = seconds_passed
        self.model_J2_effect()
        self.update_true_anomaly()


class SatelliteDynamics:
    def __init__(self, orbital_model):
        self.q = np.asarray(INITIAL_ORIENTATION, dtype=float)
        self.w = INITIAL_ANGULAR_VELOCITY.copy()
        self.I = INERTIA_MATRIX  # diagonal 3x3 matrix
        self.orbit = orbital_model

    def calculate_linear_velocity(self):
        pass
        
    def get_DCM_matrix(self):
        '''
        Create direction cosine matrix that describes transformation from ECI 
        to satellite body frame. (See page 8.)
        '''

        r = Rotation.from_quat(self.q, scalar_first=True)
        return r.as_matrix()

    def add_B_noise(self, true_B):
        static_noise = MAGNETOMETER_HARD_IRON_BIAS * 1e-6
        rng = np.random.default_rng()
        random_noise = rng.normal(loc=0.0, scale=1.8408e-6, size=3)
        least_significant_bit = 0.6e-6 # tesla per LSB
        
        unrounded_B = true_B + static_noise + random_noise
        raw_B = unrounded_B // least_significant_bit * least_significant_bit
        return raw_B

    def compute_net_torque(self, commanded_moment, B_body):
        '''
         ---
        Computes effected torque from commanded dipole moment. 
        SKIPPED DUE TO TIME CONSTRAINTS - adds disturbance torques 
        (residual atmospheric drag, gravity gradient torque, and 
        residual magnetoric dipole and Lorentz forces). Then calculates 
        total change in states.
        '''
        return np.cross(commanded_moment, B_body)

    def RK4_step_integrator(self, func, dt, t0, y0):
            '''
            The Fourth-Order Runge-Kutta (RK4) method integrates an ordinary 
            differential equation dy/dt = func(t, y), with y0 = y(t0) and
            dt = t_1 - t_0
            '''
            f1 = func(t0, y0)
            f2 = func(t0 + dt/2, y0 + (dt/2) * f1)
            f3 = func(t0 + dt/2, y0 + (dt/2) * f2)
            f4 = func(t0 + dt, y0 + dt * f3)
            y_final = y0 + (dt/6) * (f1 + 2*f2 + 2*f3 + f4)
            return y_final
    
    def update_states(self, net_torque: np.ndarray, B_rot_matrix: np.ndarray, t0: float, dt: float):
        # Euler's rigid-body equation
        def w_dot(t, y):
            return np.linalg.solve(self.I, net_torque - np.cross(y, self.I @ y))

        w_new = self.RK4_step_integrator(w_dot, dt, t0, self.w)

        # q maps ECI to body, so body angular velocity enters with a negative sign.
        w_mid = (self.w + w_new) / 2
        norm_w = np.linalg.norm(w_mid)
        inner = norm_w * dt / 2
        w1, w2, w3 = w_mid
        omega_matrix = np.array([[   0, -w1, -w2, -w3],
                                 [  w1,   0,  w3, -w2],
                                 [  w2, -w3,   0,  w1],
                                 [  w3,  w2, -w1,   0]])
        
        if norm_w == 0: # prevent division by 0
            q_new = np.asarray(self.q, dtype=float)
        else:
            q_new = (
                np.cos(inner) * np.identity(4)
                - np.sin(inner) / norm_w * omega_matrix
            ) @ self.q
            q_new /= np.linalg.norm(q_new)

        # Update state attributes
        self.w = w_new
        self.q = q_new


if __name__ == "__main__":
    main()
