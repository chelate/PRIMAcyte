using Pluto
import Pkg


pipeline = [
    "01_QualityControl.jl",
    "02_BatchCorrections.jl",
    "03_CelltypeAnnotation.jl",
    # "04_ImageViewer.jl"
    # does not generated 
    "05_tSNE.jl",
    "06_ClinicalModels.jl",
]


for file in pipeline
    println("running $file ")
    Pluto.activate_notebook_environment(joinpath(@__DIR__,file))
    include(joinpath(@__DIR__,file))
    Pkg.activate()
    println("finished $file ")
end