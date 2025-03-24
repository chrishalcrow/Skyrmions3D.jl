# Flow

After generating an initial Skyrmion, you can evolve it in a number of ways. Currently, only energy minimising flows are implemented in `Skyrmions3D`. 

## Gradient Flow

The simplest of these is a simple gradient flow, which models the differential equation

```math
\dot{\boldsymbol{\pi}} = -\frac{\delta E}{\delta \boldsymbol{\pi}}
```

Suppose we have ``B=4`` rational map skyrmion, but modified so that it doesn't have cubic symmetry (which is the symmetry of the energy minimiser).

```julia
using Skyrmions3D
skyrmion = Skyrmion(30,0.2)
p4(z) = z^4 + 3*sqrt(3)*im*z^2 + 1
q4(z) = z^4 - 3*sqrt(3)*im*z^2 + 0
make_rational_map!(skyrmion, p4, q4)
```

You can apply a gradient flow like so:

```
gradient_flow!(skyrmion)
```

When you just pass the skyrmion, the algorithm will get applied for one step with an automatically chosen time-step. If you'd like to run the gradient flow for a bit longer (advised!) or for your own custom time-step, you can as follows:

```julia
gradient_flow!(skyrmion, steps=100, dt=0.0001)

>>> initial: energy: 5.821941850733913
after 100 steps, error = 11.52
final energy: 5.545750226541662
```

Alternatively, you can specify a tolerance `tolerance`, which will allow the gradient flow to run until

```math
\text{max}_{\text{grid}} | \frac{\delta E}{\delta \boldsymbol{\pi}} | \leq \texttt{tolerance},
```

as follows:

```julia
gradient_flow!(skyrmion, tolerance=1.0)
```

Note that, depending on your boundary conditions and lattice size, your specified tolerance might be impossible to obtain.

## Arrested Newton Flow

A more efficient algorith is the Newton Flow method. This models the second order differential equation

```math
\ddot{\boldsymbol{\pi}} = -\frac{\delta E}{\delta \boldsymbol{\pi}}
```

But whenever the energy increases, the velocity is reset to zero. This has two advantages: usually, the minimum is found in less flow time; and the second-order time evolution allows for a larger time-step. Flow using Arrested Newton as follows:

```julia
arrested_newton_flow!(skyrmion, tolerance=0.01)
```
