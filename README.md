# shiny-apps
Collection of shiny apps I have created for various projects

I am happy for everyone to use these apps. Please let me know if you have used them and found them useful, and also let me know if you would like any additional features added. 

I would also appreciate it if you acknowledge me in any publications/presentations etc.

Each app will have its own subdirectory along with a README.md page with accompanying information.



# How to launch a shiny app through RStudio

### 1) Download R and RStudio
- Ensure you have R and RStudio installed on your PC. You need to have downloaded, installed, and opened R before you can install RStudio.
- Download link for R: https://cran.r-project.org/mirrors.html
- Download link for RStudio: https://posit.co/download/rstudio-desktop/

### 2) Launch RStudio

### 3) Install shiny app and other relevant packages
- You will need to install relevant packages for the app. To find out which packages these are, please see the README.md file in the apps subdirectory.
- To install the packages, go to your console and type:

```
install.packages('shiny')
install.packages('package-name1')
install.packages('package-name2')
```

### 4) Launch the app
- Once you have the shiny package and other relevant packages installed, you can now launch the app:
- You will need to change the `subdir` option to the apps subdirectory.
- In your console, type:

```
runGitHub(repo = 'Gibbatron/shiny-apps', subdir = 'app-subdirectory-name')
```
