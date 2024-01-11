# Project using Matlab

* Date : 12/2020
* Done in : 2 weeks
* Languages used : Matlab

## Purpose of the project
The program models the magnetic flux of a group of solenoids. I did some simulations to study if it is possible to focus the magnetic flux on one point so that it could improve the magnetic charging for smartphones.
I designed a partial torus using several solenoids displaced strategically. The simulation of one solenoid’s magnetic flux has been inspired by D. Cebron’s work.

Here are the parameters that can be changed:

* a: Solenoid radius

* L: Solenoid length

* Br: Residual field = field in infinite solenoid (mu*n*I)

* Theta: Angle between each solenoid

---
Here are some results I got with the simulation:

![alt text](https://github.com/Yutotto97/Focused_Magnetic_field/blob/main/Prototypemagneticflux_A%3D5.gif "video of magnetic field changing the shape")

![alt text](https://github.com/Yutotto97/Focused_Magnetic_field/blob/main/Prototype_yield_on_abscissa_A%3D5.gif "video of yield on the abscissa plotted against d")

d is a parameter that characterizes the partial torus. It is defined as shown in the picture:
![alt text](https://github.com/Yutotto97/Focused_Magnetic_field/blob/main/definition%20of%20d.png "definition of d")
