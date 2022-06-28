## Contributing

Contribute to the software

  * Setting up environment
    - Download and install [anaconda](https://www.anaconda.com/products/individual) (python version: >3.0)
    - Download and install [Qt Designer](https://build-system.fman.io/qt-designer-download) (Qt version: 5.0)

  * Preparing your own descriptor function:

    - All the core functions are arranged in separate modules. The generic structure of a module is as follows:
        
        ````python
            def your_function_name (data_stack, **vars):
                ...
                code
                ...
                return your_descriptor
                  
        ````
      data_stack : 3-D array of the polarimetric matrix (N x N x 9 (T3/C3); N x N x 4 (C2/T2)).

      \**vars :  list of required variables(E.g. **window_size**, **ellipticity** etc.)


> A template module is provided in the functions folder [mod_template](../functions/mod_template.py)

  * Updating the GUI 
    - Open the **mainWindow.ui** file in the Qt Designer and update the elements.
    - link the module with the ui elements in [SAR_Tools.py](../SAR_Tools.py)

> Please contribute by [forking](http://help.github.com/forking/) and sending a [pull request](https://docs.github.com/en/github/getting-started-with-github/github-glossary#pull).


