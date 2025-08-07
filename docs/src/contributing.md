# Contributing

There are many way to contribute to `Skyrmions3D` using the functionality of Git and GitHub. These include (but are not limited to)

* adding examples to the documentation,
* fixing bugs and typos,
* implementing new skyrmion types. 

If you want to do any of the above, first [raise an issue on the GitHub page](https://github.com/chrishalcrow/Skyrmions3D.jl/issues) to discuss the idea with others. This conversation can help understand exactly the best way to fix a bug / implement a new feature using the experience of the community. Issues should be well written: if you are reporting a bug you should provide a minimal reproducible example of the bug as well as the details of which version of the code you are using and your coding environment, if you are suggesting an enhancement provide evidence that the enhancement is worthwhile including links to relevant research. 

To then make the change, [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) the `Skyrmions3D` repository on GitHub, make the edits to your fork of the source code, test the code, and then [open a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork). After you have opened the pull request, other members of the community will review the changes. They may add comments, request tweaks to the code, or add changes of their own. More about pull requests can be read [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests). 

If you are regularly using `Skyrmions3D`, it is integral for a project you are undertaking, or the change you are planning could take a long time, you will likely want to want to use [branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches). These let you make edits to your fork of the code while continually maintaining a working version on your machine. 

Testing the code involves checking that the documentation builds and running the test suite, which your GitHub fork should do automatically (using workflows). Alternatively, you can run the tests on your local machine first. To do so, change directory to the project and [run Julia with the additional flag `--project=.`](https://pkgdocs.julialang.org/v1/environments/). If you fix a bug or add a new feature, consider adding a test to an appropriate file in the test suite (e.g. *test/do_diff_tests.jl*); more on testing can be found [in the Julia documentation](https://docs.julialang.org/en/v1/stdlib/Test/). 

When you open a pull request provide a link in the description to the existing issue this addresses, along with any other information that can better help the community evaluate the code. 

If you would like further guidance on best practices in coding, see the [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/) and [GitHub's guide for contributors](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project). 

