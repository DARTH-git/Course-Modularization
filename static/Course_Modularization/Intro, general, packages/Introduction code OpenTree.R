# Make sure you have a recent version from R
getRversion()
# You can check the newest version on the website from R
# https://www.r-project.org

# Download OpenTree
# OpenTree is currenly only available via GitHub
# To instal a package from GitHub you can use install_github
# With the repository name in it
install_github("DARTH-git/OpenTree", force = TRUE) # (Un)comment if there is a newer version

# After you are 
p_load_gh("DARTH-git/OpenTree")


## Create or open decision tree. 

#The function `create_tree()` creates a blank tree and the function `open_tree()` opens an existing decision tree. 

#*IMPORTANT*: since `create_tree()` always creates blank new tree, do not use it to access or modify an existing tree, or else the tree will get erased. Always use `open_tree()` to open and modify existing trees.

#Any changes made to the tree in OpenTree are automatically saved as a `.json` file to the working directory. If you are running it in an R script, the `.json` file will be saved to the path on your machine specified in `dir_name`. If you are running it in an R markdown document, the `.json` file will be saved to the path where the R markdown document is located. 

create_tree(file_name = "DemoTree", dir_name = getwd())
