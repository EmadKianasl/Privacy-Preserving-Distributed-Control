# Smart Building Vibration Control Using Fuzzy PID

Active vibration control of a three-story building subjected to seismic excitation using a **Fuzzy PID controller**. The project investigates structural response reduction through intelligent control, with particular emphasis on controller performance and sensor placement.

The complete model and simulations are implemented using **MATLAB/Simulink**.

---

## Overview

Buildings subjected to seismic excitation can experience significant lateral vibrations and structural responses. Active structural control provides a way to mitigate these responses by applying a control force based on measured structural states.

This project develops and evaluates an active vibration control framework for a **three-story building model** under seismic excitation.

A **Fuzzy PID controller** is implemented to regulate the structural response and reduce vibration compared with the uncontrolled building.

The project also investigates the influence of **sensor placement** on the performance of the control system.

---

## Objectives

The main objectives of this project are:

* Modeling the dynamic behavior of a three-story building
* Simulating the building response under seismic excitation
* Designing an active vibration control system
* Implementing a Fuzzy PID controller
* Investigating the effect of sensor placement
* Comparing controlled and uncontrolled structural responses
* Evaluating vibration reduction in terms of displacement and acceleration
* Developing the complete simulation framework in MATLAB/Simulink

---

## System Model

The structure is represented by a simplified **three-degree-of-freedom building model**, where each degree of freedom represents the lateral motion of a story.

The general structural dynamics can be represented in state-space form as:

$$
\dot{\mathbf{x}} =
A\mathbf{x} + B\mathbf{u} + E\mathbf{w}
$$

where:

* \(\mathbf{x}\) is the structural state vector
* \(A\) is the system matrix
* \(B\) is the control input matrix
* \(\mathbf{u}\) is the control force
* \(E\) represents the influence of the seismic excitation
* \(\mathbf{w}\) represents the external ground excitation

The controller receives structural response information through the selected sensor configuration and generates the corresponding control action.

---

## Seismic Excitation

The building is subjected to an external seismic excitation, and the structural response is evaluated both **without control** and **with active control**.

The comparison focuses on the resulting structural vibration, particularly:

* Story displacement
* Structural acceleration
* Overall vibration response

This provides a direct assessment of the effect of the active controller on the building response.

---

## Fuzzy PID Controller

The control strategy combines the conventional PID control structure with a fuzzy inference mechanism.

A conventional PID controller is based on the error:

$$
e(t)=r(t)-y(t)
$$

and its derivative:

$$
\dot e(t)=\frac{de(t)}{dt}
$$

The PID control law can be expressed as:

$$
u(t)
=
K_Pe(t)
+
K_I\int e(t)\,dt
+
K_D\frac{de(t)}{dt}
$$

In the Fuzzy PID approach, fuzzy logic is used to determine or adjust the controller behavior based on the system error and its variation.

This provides a nonlinear control mechanism that can adapt the control action according to the instantaneous structural response.

---

## Sensor Placement

An important part of the project is the investigation of **sensor placement**.

The location of the measured structural response can influence the information available to the controller and consequently affect the control performance.

Different sensor configurations are therefore considered to investigate their effect on:

* Structural displacement
* Structural acceleration
* Controller response
* Overall vibration reduction

---

## Simulation Framework

The project is implemented using:

**MATLAB**

and

**Simulink**

The simulation workflow consists of:

```text
Seismic Excitation
        │
        ▼
┌──────────────────┐
│  Building Model  │
│    3 Stories     │
└────────┬─────────┘
         │
         ▼
   Structural Response
         │
         ▼
   Sensor Measurement
         │
         ▼
┌──────────────────┐
│   Fuzzy PID      │
│    Controller    │
└────────┬─────────┘
         │
         ▼
    Control Force
         │
         └──────────────► Building
```

---

## Results

The controlled and uncontrolled responses are compared to evaluate the effectiveness of the proposed control strategy.

According to the results reported for this project, the active Fuzzy PID control strategy can achieve **up to approximately 85% reduction in displacement and acceleration responses** compared with the uncontrolled case, depending on the response quantity and sensor configuration.

The results demonstrate the potential of fuzzy PID control for reducing seismic-induced vibration in the modeled three-story structure.

> **Note:** The reported reduction represents the results of the specific simulation cases included in this project and should not be interpreted as a universal performance guarantee for other building models or seismic records.

---

## Project Structure

The current repository contains the main project implementation inside:

```text
smart control of building (Vibration Control)/
```

The directory contains the MATLAB/Simulink files associated with the building vibration-control study.

A simplified representation of the repository is:

```text
smart-control-of-building-Vibration-Control/
│
├── smart control of building (Vibration Control)/
│   ├── MATLAB files
│   ├── Simulink models
│   └── Supporting project files
│
└── README.md
```

---

## Software

The project was developed using:

* **MATLAB**
* **Simulink**
* Fuzzy Logic control modeling
* Dynamic system simulation
* Structural vibration analysis

---

## Key Topics

This project covers several topics in structural dynamics and control:

* Structural Dynamics
* Earthquake Engineering
* Seismic Response
* Vibration Control
* Active Structural Control
* Fuzzy Logic Control
* PID Control
* Fuzzy PID Control
* Sensor Placement
* State-Space Modeling
* MATLAB
* Simulink

---

## How to Run

### Requirements

* MATLAB
* Simulink
* Required MATLAB/Simulink toolboxes for the included models

### Workflow

1. Clone the repository:

```bash
git clone https://github.com/EmadKianasl/smart-control-of-building-Vibration-Control.git
```

2. Open MATLAB.

3. Navigate to:

```text
smart control of building (Vibration Control)/
```

4. Open the corresponding Simulink model.

5. Run the simulation and compare the controlled and uncontrolled structural responses.

---

## Results and Visualization

The project produces time-history responses that can be used to compare the structural behavior before and after applying active vibration control.

Typical outputs include:

* Story displacement response
* Story acceleration response
* Controlled response
* Uncontrolled response
* Controller response
* Comparison of different sensor configurations

---

## Academic Context

This repository contains an academic engineering project focused on the application of intelligent control techniques to structural vibration problems.

**Application:** Active vibration control of building structures
**Structure:** Three-story building model
**Excitation:** Seismic excitation
**Controller:** Fuzzy PID
**Software:** MATLAB / Simulink

---

## Author

**Emad Kian Asl**

Mechanical Engineering

Isfahan University of Technology

---

## License

This repository is provided primarily for **academic and educational purposes**.

Please refer to the repository files for the original MATLAB/Simulink implementation and project materials.
