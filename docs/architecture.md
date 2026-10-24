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
- **Dynamics**: articulations, and the algorithms that move them. Small machines go through their
  mass matrix; large ones through the articulated-body recursion, in time linear in their size.
- **Solver**: each machine's loops, drives, limits and springs. Pistons and servos hold like
  hydraulics: rigid once their valve shuts, and a released lever stops them where they are.
- **Contacts**: contacts, friction, and resistance to rolling and spinning, solved through the
  machines they touch. A contact reported before it touches is speculative, so a fast body stops
  on a surface instead of passing through it. A contact the host marks as part of a mechanism, a
  cam on its follower or a gear on its gear, is solved with its machine's joints instead, exactly,
  its friction bounded by however hard it ends up pushing. Each island is solved on its own
  budget: one holding a machine as hard as a machine needs, one of rigid bodies alone on lean rows
  made once a step, its impacts swept harder than its resting rows.
- **World**: what the host talks to. It takes bodies, joints and contacts, works out the
  machines, and steps them. A body is posed by its origin, its centre of mass wherever its mass
  is, so an edit that moves the mass leaves the body and its joints where they were. A free body
  added or edited changes its own group, and a joint its own machines: the world is rebuilt only
  for a weld.
- **Moves**: one call for each way the host moves a body. A machine moved whole keeps its joints;
  a static anchor moved carries what hangs on it; a kinematic body is moved by the host alone and
  pushes and carries what it touches, a machine anchored on it riding on its motion; a joint is
  set by hand from the machine's root; a grab pulls a point toward the mouse, never harder than
  its cap.
- **Sleep**: a body still for a few steps falls asleep and stops being stepped, until a push, an
  edit or an awake neighbour wakes it. A machine sleeps whole. What is stepped splits into
  islands, the bodies touching or jointed to each other. Contacts between sleepers stay in the
  world with their impulses, so the host hands only pairs with an awake end, and a sleeper woken
  stands on what held it.
- **Faults**: a row only writes into the moving bodies of its own island, never into the ground
  or a sleeper it touches, so a NaN born in an island stays there. Once the step is done the world
  looks for one: an island that has one, or a body faster than the host allows, is put back as the
  step found it, at rest, and the host is told. Bad input never gets that far: a call given a
  value the world can't use refuses it.

Heavy Feather does no collision detection. The host finds contacts with its own, hands them
over, steps the world, and reads the poses back.

## Why

- **Featherstone**, because heavy machinery needs joints that hold under load and real masses.
  In the prototype it did both, and it was faster than the alternatives on machine scenes.
- **Double precision**, because the solver relies on numbers too small for float.
- **No dependencies**, so it fits into any engine.
- **No collision detection**, because every engine already has its own, tuned to its geometry. (A general purpose collision detection might get provided later)
