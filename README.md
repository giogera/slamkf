# slamkf
Simultaneous localization and mapping (SLAM) using Kalman filtering on LIDAR data acquired from an Omron Adept MobileRobots Pioneer LX mobile robot.

This project was created for an exam, so it is not maintained. 

## Problem definition
The robot moves in a 2D environment described by a $(x, y)$ Cartesian reference frame. The robot pose at time $t$ is defined by its spatial coordinates $x(t)$, $y(t)$, and its orientation $θ(t)$ with respect to the $x$ axis. The robot motion model is described by the discrete-time equations
```math
    \begin{aligned}
        x(t + 1) &= x(t) + T_s u_f (t) \cos \theta(t) \\
        y(t + 1) &= y(t) + T_s u_f (t) \sin \theta(t) \\
        θ(t + 1) &= \theta(t) + T_s u_a(t)
    \end{aligned}
```
where $u_f(t)$ and $u_a(t)$ are respectively the forward and angular velocities at time $t$, and $T_s$ is the sampling time. In this model, we suppose the odometric measurements are corrupted by a white process $w(t)$ with zero mean and covariance matrix
```math
    Q(t) = 
    \begin{bmatrix}
        \sigma_f^2(t) & 0 \\
        0 & \sigma_a^2(t) \\
    \end{bmatrix},
```
which is supposed to be directly added to the nominal values; moreover, $Q(t)$ changes according to a threshold on the angular velocity. The robot is equipped with a Lidar scanner which collects range measurements at predefined angles in clockwise sense. At each time step $t$, the Lidar scanner gives us a range value $\rho_k(t)$ and an angle value $\alpha_k(t)$ for each scan $k$, measured with respect to the robot’s current pose. The measurements are also corrupted by an additive white noise $v(t)$ with zero mean and covariance matrix
```math
    R = 
    \begin{bmatrix}
        \sigma_\rho^2 & 0 \\
        0 & \sigma_\alpha^2 \\
    \end{bmatrix}.
```
Let $z(t)$ be a state variable containing the pose of the robot and the Cartesian coordinates of the landmarks. An Extended Kalman Filter (EKF) has been implemented to correct the pose of the robot using odometric data provided by the sensors and the measurements related to the landmarks spotted at each time step $t$.

## Repository structure
```
.
├── lidar_scan.mat              # Lidar scan data
├── README.md                   # The file you are reading right now
└── slam.m                      # Main script
```

## Usage
1. Clone the repo:
    ```
    git clone https://github.com/giogera/slamkf.git
    ```
2. Open `slam.m` in MATLAB and run it.
