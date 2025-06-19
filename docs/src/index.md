# Skyrmions3D

A [Julia](http://julialang.org) package for creating, manipulating and probing Skyrmions, in the Skyrme model of nuclear physics.

---

## Installation

If you don't have it, install [Julia](https://docs.julialang.org/en/v1/manual/installation/).

Once Julia is installed, you can enter the language REPL by typing `julia` into your terminal / command line. When inside, you should see `julia>` at the left-hand side of your terminal.

Now go into package mode by typing `]`. Now something like `@(v1.11) >` should appear at the left-hand side of your terminal. You are now in Julia's package manager. From here, you can install `Skyrmions3D` by typing

`add https://github.com/chrishalcrow/Skyrmions3D.jl.git`

When installing the package, Julia will install all other packages `Skyrmions3D` depends on. This might take a little while.

Once it's installed, go back to the Julia REPL (either by typing backspace from package mode, or by typing `julia` into your base terminal). You can now "use" Skyrmions3D by typing `using Skyrmions3D`. If this works, you've installed the package!

If you have any problems installing the package, please post an [issue on the github page](https://github.com/chrishalcrow/Skyrmions3D.jl/issues).

### Learning resources

Details of how to make, transform, compute properties of, flow and visualise skyrmions can be found in the sidebar. The API contains descriptions of all exposed functions from the package.

There is  a tutorial available on the main [github](https://github.com/chrishalcrow/Skyrmions3D.jl.git) repo, either as a Julia file or a Jupyter notepad.

I have also made [one](https://youtu.be/TI5huk6Rqos) or [two](https://youtu.be/HeSs7yVGXR4) YouTube videos about the package.

### Authors

I am Chris Halcrow, a research software engineer at the University of Edinburgh. I would love there to be more authors of this package. Please join in.

### Contributing

There are many way to contribute to `Skyrmions3D`. These include

* adding examples to the documentation,
* fixing bugs and typos,
* implementing new skyrmions types. 

If you want to do any of the above, first [raise an issue on the github page](https://github.com/chrishalcrow/Skyrmions3D.jl/issues) to discuss the idea with others. To then make the change, fork the `Skyrmions3D` repository on github, make the edits to your fork of the code, and then open a pull request. Lots of details of how to complete this process are available [from GitHub](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork). 

 
