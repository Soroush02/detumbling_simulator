## Announcement: The code is (partially) functional!

See it in action below:

This is a 1-hour simulation:
<img width="640" height="480" alt="Figure_1" src="https://github.com/user-attachments/assets/92373a65-f2d0-4c2e-80f1-f4e4296a544a" />

A full simulation yields unexpected results:
<img width="640" height="480" alt="Figure_2" src="https://github.com/user-attachments/assets/27749ca3-d003-4a02-82be-97c2aee5d5e6" />
(A full 8-hour simulation takes around 1 hour to compute.)

#### Issues fixed the morning after submission:

* **Changed starting parameters** - I realised I entered random placeholders for some parameters (e.g. inertia tensors = 1) and forgot to change them to realistic values.
* Changed line 155 from `np.min(INERTIA_MATRIX)` to `np.min(INERTIA_MATRIX.diagonal())`. Since the former always output 0 (for a diagonal matrix), `BDot.update_desired_moment()` would always output `[0. 0. 0.]`; this is now fixed.
* Deleted method `SimulationController.transform_B_Body()`, which repeatetd code already in `SimulationController.stepchange()`.
* Make progress indicator show every 5% for clarity.
* Add `plt.savefig()` (line 122) to save simulation results.
