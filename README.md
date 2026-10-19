## Announcement: The code is fully functional!

See it in action below! (This is only a 1-hour simulation; a full 8-hour simulation takes around 1 hour to compute.)
<img width="640" height="480" alt="Figure_1" src="https://github.com/user-attachments/assets/92373a65-f2d0-4c2e-80f1-f4e4296a544a" />

#### Issues fixed the morning after submission:

* Change starting parameters - I realised I entered random placeholders for some parameters (e.g. inertia tensors 1) and forgot to change them to realistic values.
* Changed line 155 from `np.min(INERTIA_MATRIX)` to `np.min(INERTIA_MATRIX.diagonal())`. Since the former always output 0 (for a diagonal matrix), BDot.update_desired_moment() would always output [0. 0. 0.]; this is now fixed.
* Deleted method SimulationController.transform_B_Body(), which duplicated code in SimulationController.stepchange().
* Make progress indicator show every 5% for clarity.
