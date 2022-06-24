<div style="font-family:Courier">

# FROM THE OUTPUT OF THE RUN_MOTORS CODE
***2022 06 23***

## From the "out" variable:
- **Cols 1 to 6 = pitching and heaving positions for the three rigs (starting with the foremost upstream)**
    01: Shawn   pitching (defunct)
    02: Shawn   heaving  (defunct)
    03: Gromit  pitching (leading) [rad]
    04: Gromit  heaving  (leading) [m]
    05: Wallace pitching (trailing) [rad]
    06: Wallace heaving  (trailing) [m]
- **Cols 13 to 16 = vectrino data.**
    13: x flow velocity [m/s]
    14: y flow velocity [m/s]
    15: z1 flow velocity [m/s]
    16: z2 flow velocity [m/s]
- **Cols 17 to 22 = Leading foil data:**
    17: normal force (Fy) [N]
    18: tangential force (Fx) [N]
    22: spanwise torque (Tz) [N]
- **Cols 7 to 12 = Trailing foil data:**
    07: normal force (Fy) [N]
    08: tangential force (Fx) [N]
    12: spanwise torque (Tz) [N]

## From the "Prof_out_angle" variable:
- **Cols 1 to 6 = pitching and heaving positions for the three rigs (starting with the foremost upstream)**
    01: Shawn   pitching (defunct)
    02: Shawn   heaving  (defunct)
    03: Gromit  pitching (leading) [deg]
    04: Gromit  heaving  (leading) [m]
    05: Wallace pitching (trailing) [deg]
    06: Wallace heaving  (trailing) [m]

</body>