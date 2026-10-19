## Update: The Code is Fully Functional!

See it in action below:

<p align="center">
  <img width="auto" height="280px" alt="Figure_gg" src="https://github.com/user-attachments/assets/7f2cd89f-b433-434b-a49f-5ecb971e76bd" />
  <img width="auto" height="280px" alt="Figure_Final" src="https://github.com/user-attachments/assets/11d7fbf9-e804-433c-ba6d-1bef483d35fc" />
</p>

Time taken to compute scales proportionally to `MAX_SIM_HOURS` and inversely proportionally to `DT_SIMULATOR`. For reference, a full 8-hour simulation takes around 1 hour to compute with `DT_SIMULATOR = 0.05`.

#### Issues fixed the day after submission

* **Changed starting parameters** - I realised I had entered random placeholders for some parameters (e.g. inertia tensors = 1) and forgot to update them to realistic values.
* Changed the minus in line 475 (`qnew = (... - ...)`) to a plus, fixing a typo to make it consistent with my logbook.
* Changed line 155 from `np.min(INERTIA_MATRIX)` to `np.min(INERTIA_MATRIX.diagonal())`. Since the former always outputs 0 (for a diagonal matrix), `BDot.update_desired_moment()` would always output `[0. 0. 0.]`; this is now fixed.
* Added the `.inv()` modifier to the rotation defined in `SatelliteDynamics.get_DCM_matrix()`. This ensures the rotation is passive rather than active.
* Added a check in `BDotController.stepchange()` to prevent unnecessary commanded moments if angular velocity is already below 1º/sec.
* Deleted method `SimulationController.transform_B_Body()`, which repeated code already in `SimulationController.stepchange()`.
* Make the progress indicator show every 1% for clarity.
* Add `plt.savefig()` (line 122) to save simulation results.
