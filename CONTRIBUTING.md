# Contributing to PROBAB3
The open source license for this repository is the [MIT license](LICENSE.md).

Bug fixes and new features are more than welcome. 

## Contribution Process

To insure an orderly and efficent contribution process, please follow the following outlined steps.

1. Install a development enviornment as instructed in the [README](README.md#development-enviornment-installation)

2. Modify the code as you please locally.

3. Test backwards compatibility

    After changing or adding a feature to the code, you must make sure other functionality did not break. 

    1. Run the existing repo tests with pytest 
        ```
        poetry run pytest -vvv
        ```
        
        Make sure they all passed.
        Save a picture of this for later.
    
    2. Manually run the whole [check_sampling_dists](./probab3/tests/check_sampling_dists.ipynb) notebook.

4. Add tests to your feature

    To ensure maintainability of features and bug fixes, please add a simple unittest for your feature. This will make sure that future modifications will not break your feature by accident.

5. Commit your changes to a new [feature branch](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging)

6. Open a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) on a seperate feature branch 

    Make sure your pull request has the following:
    
    1. A desciption of your feature
    
    2. A picture of the passed tests on your local computer

    3. New run of the check_sampling_dists file (should be committed along with your changes)

    4. A line containing a description of the addition in the CHANGELOG.md file

7. Add the owner of the repository as a reviewer for your code. Contact her by email to let her know you are waitning on a code review.

8. After an approval, you may be able to merge your changes to the main branch.

9. You would be able to see your changes in the following version release after your merger.

