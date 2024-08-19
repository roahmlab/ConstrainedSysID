# System Identification For Constrained Robots

## Introduction

Our paper is available [here](https://arxiv.org/abs/2408.08830).

Identifying the parameters of robotic systems, such as motor inertia or joint friction, is critical to satisfactory controller synthesis, model analysis, and observer design. 
Conventional identification techniques are designed primarily for unconstrained systems, such as robotic manipulators. 
In contrast, the growing importance of legged robots that feature closed kinematic chains or other constraints, poses challenges to these traditional methods. 
This paper introduces a system identification approach for constrained systems that relies on iterative least squares to identify motor inertia and joint friction parameters from data.
The proposed approach is validated in simulation and in the real-world on Digit, which is a 20 degree-of-freedom humanoid robot built by Agility Robotics.
In these experiments, the parameters identified by the proposed method enable a model-based controller to achieve better tracking performance than when it uses the default parameters provided by the manufacturer. 

![Summary Figure](https://github.com/user-attachments/assets/90796e01-97d6-4b8c-89a1-65b0e4341259)

To illustrate the utility of this method, this paper compared the performance of a model based tracking controller to track a user-specified trajectory. 
The tracking performance of the controller was significantly less when using the parameters identified by the algorithm developed in this paper (drawn in red on the right) when compared to using the parameters specified by the manufacturer (drawn in blue on the right).
The actual behavior of the robot while following the user-specified trajectory using the parameters identified by the algorithm developed in this paper can be seen on the top right of the image at 4 time instances. 

## Credits

Daniel Haugk (st161112@stud.uni-stuttgart.de)

Bohao Zhang (jimzhang@umich.edu)

This work is developed under [RoahmLab](https://www.roahmlab.com/), University of Michigan, Ann Arbor.

To cite our work in your academic research, please use the following bibtex entry:

```bibtex
@misc{zhang2024identificationconstrainedrobots,
      title={System Identification For Constrained Robots}, 
      author={Bohao Zhang and Daniel Haugk and Ram Vasudevan},
      year={2024},
      eprint={2408.08830},
      archivePrefix={arXiv},
      primaryClass={cs.RO},
      url={https://arxiv.org/abs/2408.08830}, 
}
```
