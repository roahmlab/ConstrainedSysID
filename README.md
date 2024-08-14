# ConstriainedSysID: System Identification For Constrained Robots

Identifying the parameters of robotic systems, such as motor inertia or joint friction, is critical to satisfactory controller synthesis, model analysis, and observer design. 
Conventional identification techniques are designed primarily for unconstrained systems, such as robotic manipulators. 
In contrast, the growing importance of legged robots that feature closed kinematic chains or other constraints, poses challenges to these traditional methods. 
This paper introduces a system identification approach for constrained systems that relies on iterative least squares to identify motor inertia and joint friction parameters from data.
% Despite sophisticated system identification approaches, there could still be a possible mismatch between the identified parameters and the ground truth, which downgrades the tracking performance of a model-based controller.
The proposed approach is validated in simulation and in the real-world on Digit, which is a 20 degree-of-freedom humanoid robot built by Agility Robotics.
In these experiments, the parameters identified by the proposed method enable a model-based controller to achieve better tracking performance than when it uses the default parameters provided by the manufacturer. 

# Credits

Daniel Haugk (st161112@stud.uni-stuttgart.de)
Bohao Zhang (jimzhang@umich.edu)

[RoahmLab](https://www.roahmlab.com/), University of Michigan, Ann Arbor
