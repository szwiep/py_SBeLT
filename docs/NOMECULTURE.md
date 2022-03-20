# Nomeculture

## A quick overview of terms used in source code and documentation

- **Stream**: An 2D domain that consists of bed particles and model particles and is divided by subregions. The left-most point of the stream is considered the _upstream boundary_ while the right-most is considered the _downstream boundary_.
- **Bed Particle**: A particle comprising the stream bed. Bed particles do not move, they are static and fixed in their position in the bed. They represent the lowest elevation in the stream.
- **Model Particle**: A particle that can move, i.e be entrained. When not entrained, model particles rest on top of bed particles. Model particles may stack on top of other model particles. When a model particle exceeds the stream's downstream boundary it is immediately sent to the downstream boundary.
- **Supporting Particle**: A particle that directly 'holds up' model particles when they stack. Both bed and model particles can be considered supporting particles.
- **Event Particle**: A particle that has been selected for entrainment. Only model particles can be considered event particles.
- **Available Vertex**: A horizontal location in the stream where a model particle is permitted to be deposited. Available vertices are created by two particles of the same elevatino touching. See Figure 1 in [PAPER.md](https://github.com/szwiep/py_SBeLT/blob/master/paper/paper.md) for visual intuition.
- **Desired Hop**: A location in the stream that a particle selected for entrainment is set to move to based on sampling from a log-Gaussian or Gaussian (see [THEORY.md](https://github.com/szwiep/py_SBeLT/blob/master/THEORY.md)). A desired hop is not always equal to the movement implemented by the model because model particles can only be deposited at available vertices.