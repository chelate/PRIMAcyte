using Pluto

# List of notebooks to be run in order
# notebooks = filter(readdir(@__DIR__)) do x
#     num = split(x,"_")[1]
#     in(["01","02","03","04"])(num)
# end

# # Function to open and run each notebook
# for notebook_path in notebooks
#     Pluto.run(notebook = notebook_path)
# end

Pluto.run()
