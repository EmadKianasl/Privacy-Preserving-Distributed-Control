# Privacy-Preserving Distributed Control

A privacy-preserving distributed control framework for a **2-DOF vibrating mechanical system** using an **encrypted PID controller** and **Shamir's Secret Sharing** for secure multiparty computation.

The project investigates how control computation can be performed in a distributed/cloud-based environment while reducing the exposure of sensitive system information.

Implementations and supporting codes are provided in **Python/Jupyter Notebook, SageMath, and MATLAB**.

---

## Overview

Modern control systems increasingly rely on distributed and cloud-based computation. While this architecture provides computational and communication advantages, transferring system states and control-related information to external computing nodes can introduce privacy and security concerns.

This project studies a privacy-preserving control architecture in which sensitive control information is protected using cryptographic techniques while the controller operates in a distributed computational environment.

The proposed framework combines:

* Distributed control
* Encrypted control computation
* Shamir's Secret Sharing
* Secure Multiparty Computation (SMPC)
* PID control
* Numerical simulation of mechanical vibration
* Newmark time-integration method

The framework is demonstrated on a **two-degree-of-freedom vibrating system**.

---

## Key Features

* 🔐 **Privacy-preserving control computation**
* 🤝 **Secure multiparty computation**
* 🔑 **Shamir's Secret Sharing**
* 🎛️ **Encrypted PID controller**
* ⚙️ **2-DOF mechanical vibration model**
* 📐 **Newmark numerical integration**
* 🐍 **Python / Jupyter implementation**
* 🔢 **SageMath cryptographic computation**
* 🧮 **MATLAB dynamic-system simulation**

---

## System Architecture

The overall concept can be summarized as:

```text
             Physical System
                   │
                   │ System States
                   ▼
          ┌─────────────────┐
          │   Measurement   │
          └────────┬────────┘
                   │
                   ▼
          ┌─────────────────┐
          │ Privacy /       │
          │ Secret Sharing  │
          └────────┬────────┘
                   │
        ┌──────────┼──────────┐
        │          │          │
        ▼          ▼          ▼
     Node 1     Node 2     Node N
        │          │          │
        └──────────┼──────────┘
                   │
                   ▼
       Secure Distributed
        Control Computation
                   │
                   ▼
          Encrypted PID
             Controller
                   │
                   ▼
             Control Input
                   │
                   ▼
          ┌─────────────────┐
          │  2-DOF System   │
          │   Dynamics      │
          └─────────────────┘
```

The main objective is to perform the required control computation without directly exposing the sensitive system information to individual computational parties.

---

## Mechanical System

The controlled plant is modeled as a **two-degree-of-freedom vibrating mechanical system**.

The generalized coordinates can be represented as:

$$
\mathbf{x}
=
\begin{bmatrix}
x_1 \\
x_2
\end{bmatrix}
$$

with corresponding velocities:

$$
\dot{\mathbf{x}}
=
\begin{bmatrix}
\dot{x}_1 \\
\dot{x}_2
\end{bmatrix}
$$

The system dynamics can generally be represented in matrix form as:

$$
M\ddot{x}+C\dot{x}+Kx=Bu
$$

where:

* \(M\) is the mass matrix
* \(C\) is the damping matrix
* \(K\) is the stiffness matrix
* \(x\) is the displacement vector
* \(u\) is the control input
* \(B\) is the control-input matrix

The numerical response of the system is evaluated using the **Newmark method**.

---

## PID Control

The baseline controller is based on the classical PID structure.

For a tracking error \(e(t)\), the control input is expressed as:

$$
u(t)
=
K_Pe(t)
+
K_I\int e(t)\,dt
+
K_D\frac{de(t)}{dt}
$$

where:

* \(K_P\) is the proportional gain
* \(K_I\) is the integral gain
* \(K_D\) is the derivative gain

In this project, the control computation is combined with privacy-preserving cryptographic operations.

---

## Privacy-Preserving Control

A central component of the project is the use of **Shamir's Secret Sharing**.

Instead of directly transmitting sensitive control variables to a single computational entity, secret information can be divided into multiple shares.

Conceptually:

```text
              Secret Value
                   │
                   ▼
          ┌─────────────────┐
          │ Secret Sharing  │
          └────────┬────────┘
                   │
        ┌──────────┼──────────┐
        ▼          ▼          ▼
      Share 1    Share 2    Share N
        │          │          │
        ▼          ▼          ▼
     Party 1    Party 2    Party N
        │          │          │
        └──────────┼──────────┘
                   │
                   ▼
          Secure Computation
                   │
                   ▼
          Reconstructed Result
```

The purpose is to allow distributed computation while preventing any individual party from obtaining the original secret directly.

---

## Secure Multiparty Computation

The cryptographic layer is designed around the concept of **Secure Multiparty Computation (SMPC)**.

The participating computational nodes jointly perform the required operations while keeping the underlying sensitive information protected.

The project therefore connects two traditionally separate areas:

**Control Engineering**

$$
\text{System Modeling}
\rightarrow
\text{State Measurement}
\rightarrow
\text{PID Control}
$$

and

**Privacy-Preserving Computation**

$$
\text{Secret Sharing}
\rightarrow
\text{Distributed Computation}
\rightarrow
\text{Secure Control Result}
$$

The combination provides a framework for investigating privacy-aware distributed control systems.

---

## Numerical Integration

The mechanical dynamics are simulated using the **Newmark method**, a widely used numerical time-integration approach for structural and mechanical dynamic systems.

The numerical simulation is used to evaluate the time-domain response of the 2-DOF system under the applied control strategy.

---

## Software and Tools

The repository contains implementations across multiple computational environments:

| Environment          | Purpose                                          |
| -------------------- | ------------------------------------------------ |
| **MATLAB**           | Dynamic-system modeling and numerical simulation |
| **Python / Jupyter** | Computational experiments and analysis           |
| **SageMath**         | Cryptographic and secret-sharing computations    |

---

## Repository Structure

```text
Privacy-Preserving-Distributed-Control/
│
├── Privacy-Preserving-Distributed-Control/
│   ├── Python / Jupyter notebooks
│   ├── SageMath implementations
│   ├── MATLAB scripts
│   └── Supporting files
│
└── README.md
```

The implementations are organized around the different computational components of the privacy-preserving control framework.

---

## Workflow

A typical computational workflow is:

```text
1. Define the mechanical system
             │
             ▼
2. Model the 2-DOF dynamics
             │
             ▼
3. Define the PID controller
             │
             ▼
4. Protect sensitive variables
   using secret sharing
             │
             ▼
5. Perform distributed /
   secure computation
             │
             ▼
6. Obtain the control input
             │
             ▼
7. Integrate system dynamics
   using Newmark method
             │
             ▼
8. Evaluate system response
```

---

## Research Topics

This project lies at the intersection of several research areas:

* **Distributed Control**
* **Privacy-Preserving Control**
* **Cyber-Physical Systems**
* **Secure Multiparty Computation**
* **Secret Sharing**
* **Control Systems**
* **PID Control**
* **Mechanical Vibrations**
* **Structural Dynamics**
* **Cloud-Based Control**
* **Cryptographic Computing**
* **Numerical Time Integration**

---

## Motivation

Cloud-based and distributed control architectures can provide computational flexibility and scalability, but they may also require sensitive system information to be shared across computational nodes.

Privacy-preserving techniques provide a potential approach for reducing direct exposure of such information.

This project explores this concept through the integration of:

$$
\boxed{
\text{Distributed Control}
+
\text{Cryptography}
+
\text{Mechanical Dynamics}
}
$$

with a 2-DOF vibrating system serving as the computational case study.

---

## Limitations

This repository represents a computational study of privacy-preserving distributed control using a simplified mechanical system.

The results should therefore be interpreted within the assumptions and modeling choices implemented in the provided MATLAB, Python/Jupyter, and SageMath codes.

The framework is intended as a research and educational implementation rather than a production-ready secure control platform.

---

## Getting Started

### Requirements

Depending on the part of the project being used, the following software may be required:

* MATLAB
* Python
* Jupyter Notebook
* SageMath

Additional MATLAB toolboxes or Python packages may be required by individual scripts or notebooks.

### Clone the Repository

```bash
git clone https://github.com/EmadKianasl/Privacy-Preserving-Distributed-Control.git
```

Then navigate to the project directory:

```bash
cd Privacy-Preserving-Distributed-Control
```

Open the relevant MATLAB scripts, Jupyter notebooks, or SageMath files according to the component you want to reproduce.

---

## Academic Scope

This project is focused on the integration of **control engineering and privacy-preserving computation**.

It can serve as a starting point for further investigation of:

* Privacy-preserving feedback control
* Secure cloud-based control
* Distributed controller architectures
* Cryptographic control systems
* Secure cyber-physical systems
* Privacy-aware multi-agent control

---

## Author

**Emad Kian Asl**

Mechanical Engineering

Isfahan University of Technology

---

## License

This repository is provided primarily for **academic and educational purposes**.

Please refer to the source files for the implementation details and computational experiments.
