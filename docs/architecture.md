# Architecture

Heavy Feather grew out of the physics prototype in Heavy Metal. This is a rough picture of how
it fits together.

## The idea

Bodies connected by joints form a **machine**, and each machine is simulated as one Featherstone
articulation. Its joints can't drift apart, because each part's pose comes from its parent's
instead of being pulled back toward it. Whatever the tree can't express, like closed loops,
motors, limits and springs, is solved directly for the whole machine at once.

Contacts are solved against the whole machine too, so a wheel touching the ground feels the
mass of everything attached to it.

## The systems

- **Math**: vectors, quaternions and spatial algebra, all in double precision.
- **Dynamics**: articulations, and the algorithms that move them.
- **Solver**: each machine's loops, drives, limits and springs.
- **Contacts**: contacts and friction, solved through the machines they touch.
- **World**: what the host talks to. It takes bodies, joints and contacts, works out the
  machines, and steps them.

Heavy Feather does no collision detection. The host finds contacts with its own, hands them
over, steps the world, and reads the poses back.

## Why

- **Featherstone**, because heavy machinery needs joints that hold under load and real masses.
  In the prototype it did both, and it was faster than the alternatives on machine scenes.
- **Double precision**, because the solver relies on numbers too small for float.
- **No dependencies**, so it fits into any engine.
- **No collision detection**, because every engine already has its own, tuned to its geometry. (A general purpose collision detection might get provided later)
