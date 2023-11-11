# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   usethis::use_test("name")
#   Install chocolate using windows powershell
#   Rcmd choco install git.install
#   usethis::use_git_config(user.name = "Yidi Deng", user.email = "meiosis97@gmail.com")
#   usethis::create_github_token()
#   git config --global --list
#   git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY.git (Generate a colon to the local machine)
#   cd myrepo (Treat the github repository as a directory)
#   echo "A line I wrote on my local computer  " >> README.md (An example modification of the repository)
#   git status (Check what modifications has been made)
#   git add README.md (Add modifications that hopes to be committed on the local machine.)
#   git commit -m "A commit from my local computer" (Commit the added changes to the local machine)
#   git push (Synchronize the local modification to the github repository)

hello <- function() {
  print("Hello, world!")
}
