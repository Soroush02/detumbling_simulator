## Update: The code is (somewhat) functional

See it in action below:

This is a 2-hour simulation (with `DT_SIMULATOR = 0.02`):

<img width="640" height="480" alt="Figure_3" src="https://github.com/user-attachments/assets/df8aba1a-b61c-4a0b-a77d-bcd8f827e726" />

A full simulation yields unexpected results:

<img width="640" height="480" alt="Figure_2" src="https://github.com/user-attachments/assets/27749ca3-d003-4a02-82be-97c2aee5d5e6" />

(This full 8-hour simulation takes around 1 hour to compute with `DT_SIMULATOR = 0.05`)

#### Issues fixed the morning after submission:

* **Changed starting parameters** - I realised I entered random placeholders for some parameters (e.g. inertia tensors = 1) and forgot to change them to realistic values.
* Changed line 155 from `np.min(INERTIA_MATRIX)` to `np.min(INERTIA_MATRIX.diagonal())`. Since the former always output 0 (for a diagonal matrix), `BDot.update_desired_moment()` would always output `[0. 0. 0.]`; this is now fixed.
* Deleted method `SimulationController.transform_B_Body()`, which repeatetd code already in `SimulationController.stepchange()`.
* Make progress indicator show every 5% for clarity.
* Add `plt.savefig()` (line 122) to save simulation results.
