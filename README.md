Two numerical experiments demonstrate the passive/active fault identification algorithms.


## Mars Satellite Experiment

[View the PDF for detailed setup of the experiment](test_mars_sat/doc/Satellite_Simulation_Setup.pdf)

This system has 12 states: the 3D cartesian position and velocity, the three-coordinate attitude defined by Modified Rodrigues Parameters (MRP), and the angular rate. Because the control inputs are the torque along the satellite’s principal axes, only the MRP and angular rate are used for fault identification, while position and velocity evolve with time. However, the attitude and angular rate still depend on the satellite’s position because the pointing direction depends on its location. The measurement consists of MRP and angular rate. The main dynamics are in the +satsym package. Due to MRP’s switching law, which creates a discontinuity in the state, an unscented Kalman filter (UKF) is used for estimation. 

[![Watch the video](https://raw.githubusercontent.com/yourusername/yourrepository/main/assets/thumbnail.jpg)](https://github.com/jordan787878/fault_id/tree/main/test_mars_sat/figs/anim_two_sat.mp4)

## Two-water Tank Experiment

This system has two states: the water level of the first and two tanks. The control input (a scalar) is the inflow to the first tank. The fault hypotheses are defined by three leak rates: C1L is the leak from the first tank to outside; C12 is the leak from the first tank to the second tank; and C2L is the leak from the second tank to outside. See the previous update document for a detailed description of this system. The measurement consists of water levels of both tanks. For this simulation, we consider 5 possible fault hypotheses: 0 denotes the nominal, 1~4 denotes four distinct faults, and -1 denotes an unknown fault that has not been modeled.
