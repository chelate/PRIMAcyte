import Pluto
import Pkg

function wrap_include(path)
    println(" - evaluating $path")
    m = Module()
    Base.include(m, path)
    println(" - finished $path")
    return nothing
end


pipeline_filelist = [
    "01_QualityControl.jl", #1
    "02_BatchCorrections.jl", #2
    "03_CelltypeAnnotation.jl" #3
    #"04_ImageViewer.jl", # 4 does not generate results, is not run
    #"05_tSNE.jl", #5
    #"06_ClinicalModels.jl" #6
]

Pluto.activate_notebook_environment(joinpath(@__DIR__,pipeline_filelist[1]))
module test1
    include(joinpath(@__DIR__,"01_QualityControl.jl"))
end

Pluto.activate_notebook_environment(joinpath(@__DIR__,pipeline_filelist[2]))
module test2
    include(joinpath(@__DIR__,"02_BatchCorrections.jl"))
end

Pluto.activate_notebook_environment(joinpath(@__DIR__,pipeline_filelist[3]))
m = Module()
Base.include(m, joinpath(@__DIR__,pipeline_filelist[3]))