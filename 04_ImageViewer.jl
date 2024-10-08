### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a112be8c-2245-11ef-2c0f-1d00ed8b2f19
begin
	using PlutoUI
	using Base64
	using FileIO
	using ImageCore
	using ImageIO
	using TiffImages
	using Luxor
	using CSV
	using DataFrames
	using Memoization
	using JLD2
	import HypertextLiteral: @htl
	import Statistics: quantile
	import LinearAlgebra
	import CairoMakie
end

# ╔═╡ cbe74fb5-972f-4127-a6c7-0c50cf8522cb
function Base.show(io::IO, mime::MIME"image/png",  img::AbstractMatrix{C}) where {C<:Colorant}
        FileIO.save(FileIO.Stream{format"PNG"}(io), img)
end

# ╔═╡ 28214a08-028a-415c-b294-8a290931e8eb
begin
	struct FullImage{C<:Colorant}
		img::AbstractMatrix{C}
	end
	function Base.show(io::IO, mime::MIME"image/png",  fi::FullImage)
    	    FileIO.save(FileIO.Stream{format"PNG"}(io), fi.img)
	end
end

# ╔═╡ 9a87684f-a7e6-4301-af1c-06166546f344
html"""<style>
main {
    max-width: 70%;
    margin-left: 2%;
    margin-right: 18% !important;
}

</style>
"""

# ╔═╡ 9fcdd12a-cd55-47b4-a98a-6d1c84f4b436
TableOfContents()

# ╔═╡ 7bce7789-58d5-49e7-b0cc-f38f2829f0db
begin
	JLD2.@load joinpath(@__DIR__,"data","cluster_labels.jld2") labels cols df_clustered
	delete!(cols,"Bright")
	delete!(cols,"Dark")
	md"""
	# Load cell information
	
		- labels (cluster_index => string_name)
		- cols (string_name => RGB)
		- df (filtered data frame, does not include all segmented objects))

	logic is that the mask tiff value will be zero for the unsegmented regions (i.e. balack)
	
	`ifelse(i == 0,
			RGB(0.0),
			get(cell_coloring,i, unknown_color))`
	"""
end

# ╔═╡ 73279267-ba04-47da-9bba-81a1bc26269c
md"""
# Set marker colors
"""

# ╔═╡ 9e355daf-ada9-4519-b186-2512a63a6747
@bind color_update Button("update colors")

# ╔═╡ 714ed9c9-7307-482c-b9dc-5c5027aae3e0
md"""
# Saving / loading colors
"""

# ╔═╡ f501dad3-6a8f-43a6-9931-75cdc7451769
@bind update_save Button("load / delete color scheme")

# ╔═╡ ffd766c6-e030-484b-a769-0bd3bad19fa4
md"""
activate the definition cell above to load the scheme.
"""

# ╔═╡ 50892ac6-3480-4f1b-9f9b-bf60753989c4
@bind save_button Button("save color scheme")

# ╔═╡ a79670ce-9a64-41a8-939f-4cb4e6a84dfe
md"""
## Set patient / ROI
"""

# ╔═╡ c0170387-a1d1-4f17-a41c-ae3875c2034a
@bind cellbits MultiCheckBox(["cells","save image"])

# ╔═╡ 0aecf292-b442-4b9a-b690-06e061a6af50
@bind value_scale Slider(range(1,0,20))

# ╔═╡ ba317e1a-c67c-4f5a-b62a-56d9c0f0df86
md"""
# Zoom Box
"""

# ╔═╡ ce169cd6-91ab-4efb-9915-88d09f36038c
@bind maskbits MultiCheckBox(["masks"])

# ╔═╡ 0c29afcf-c687-4776-b870-84583727592c
md"""
# Highlight Clusters
"""

# ╔═╡ da0cc366-5471-4f93-bc6d-74739910cffb
@bind highlighted MultiCheckBox(sort(unique(df_clustered.cluster_index)))

# ╔═╡ 7d73f953-4365-433d-b7f4-0716c79805f0


# ╔═╡ b558f09a-9468-4b49-915b-a614ee009618
md"""
# _

# The code
"""

# ╔═╡ f799657d-681e-42a2-9261-f1ab1c11421f
begin 
	highlighted_counts = Dict(first(df_image.Patient) => sum(in(ci, highlighted) for ci in df_image.cluster_index) for df_image in groupby(df_clustered,:Image) )
	total_counts = Dict(first(df_image.Patient) => sum( 1 for ci in df_image.cluster_index) for df_image in groupby(df_clustered,:Image) )
end

# ╔═╡ 4d9deb1e-77e6-44ca-9d7d-fbb7101e4f39
md"""
### Get the targets from the panel
"""

# ╔═╡ 718eb743-7b7c-46dd-890b-18a51240db2e
get_bounds(x,y,width,height) = 
	(range(extrema(round.(Int, [y, y+height]))...),range(extrema(round.(Int,[x, x+width]))...))

# ╔═╡ 86841954-e590-45a9-936b-7bd653109f4c
patient_image_dict = Dict(ii[1,"Patient"]=> ii[1,"Image"] for ii in groupby(df_clustered,:Image))

# ╔═╡ 39f58298-320a-4e7e-a170-7516d85149d4
@bind selected_patient  Select(sort(collect(keys(patient_image_dict)))) #Select((first∘splitext).(image_order)) 

# ╔═╡ 08358c11-cf3f-49a1-bde9-f143899c1557
patient_img_dict = Dict(ii[1,"img_number"] => ii[1,"Patient"] for ii in groupby(df_clustered,:Image))

# ╔═╡ 3daf7f1c-4687-468d-80f8-03ede23b39e8
function encode_img_matrix(img_matrix::Matrix{RGB{Float64}})
    img = colorview(RGB, img_matrix)
    buf = IOBuffer()
    save(Stream(format"PNG", buf), img)
    seek(buf, 0)  # Reset buffer position to the beginning
    return base64encode(read(buf))
end

# ╔═╡ 952f5688-ba41-49ec-b5d7-f3e04cf9b4bc
begin
struct RectangleDrawing
    img_matrix::Matrix{RGB{Float64}}
end

	function Base.show(io::IO, m::MIME"text/html", rd::RectangleDrawing)
		img_base64 = encode_img_matrix(rd.img_matrix)
		Base.show(io, m, @htl(
		"""
		<div>
		<canvas style="position: relative; width: 100%;"></canvas>
		
		<script>
		const div = currentScript.parentElement
		const canvas = div.querySelector("canvas")
		const ctx = canvas.getContext("2d")

		function resizeCanvas() {
			canvas.width = window.innerWidth
			canvas.height = canvas.width * (img.height / img.width) // Maintain image aspect ratio
		}

		const img = new Image()
		const imgBase64 = $img_base64;
		img.src = "data:image/png;base64," + imgBase64;
		img.onload = () => {
			resizeCanvas()
			ctx.drawImage(img, 0, 0, canvas.width, canvas.height)
		}
		
		var startX, startY

		function onmove(e){
			const rect = canvas.getBoundingClientRect()
			const scaleX = canvas.width / rect.width
			const scaleY = canvas.height / rect.height
			const x = (e.clientX - rect.left) * scaleX
			const y = (e.clientY - rect.top) * scaleY

			div.value = [startX * img.width / canvas.width, startY * img.height / canvas.height, 
					(x - startX) * img.width / canvas.width , (y - startY) * img.height / canvas.height]
			div.dispatchEvent(new CustomEvent("input"))
			
			ctx.drawImage(img, 0, 0, canvas.width, canvas.height) // Redraw image to clear previous rectangle
			ctx.fillStyle = 'rgba(255, 0, 0, 0.5)' // Use transparent red for the rectangle
			ctx.fillRect(startX, startY, x - startX, y - startY)
		}
		
		canvas.onpointerdown = e => {
			const rect = canvas.getBoundingClientRect()
			const scaleX = canvas.width / rect.width
			const scaleY = canvas.height / rect.height
			startX = (e.clientX - rect.left) * scaleX
			startY = (e.clientY - rect.top) * scaleY
			canvas.onpointermove = onmove
		}
		
		canvas.onpointerup = e => {
			canvas.onpointermove = null
		}

		window.onresize = () => {
			resizeCanvas()
			ctx.drawImage(img, 0, 0, canvas.width, canvas.height)
		}
		</script>
		</div>
		"""))
	end
	md"""
	[collapsed] Code defines a javascript interface for selecting regions via the function `RectangleDrawing(image::Matrix{RGB})`
	"""
end

# ╔═╡ 6eac3e81-d949-4a59-ae77-0c16a7f29aaa
md"""
### unknown color and highlighted clusters
"""

# ╔═╡ 339eb3c7-40de-4a5d-82b2-c571cf421288
begin 
const unknown_color = RGB(.3)
const highlight = RGB(0.4,0.1,0.1)
end

# ╔═╡ 8272ba6b-b303-43ad-8b3f-a96e28e387e0
begin 
	labs = sort(collect(keys(cols)))
	push!(labs,"Unknown")
	swatches = [CairoMakie.PolyElement(color = get(cols,l,unknown_color), strokewidth = 0) for l in labs]
	cell_type_legend = CairoMakie.Figure(size = (500,250))
	CairoMakie.Legend(cell_type_legend[1, 1], swatches, labs, patchsize = (35, 35), rowgap = 10, nbanks = 3, framevisible = false)
	cell_type_legend
end

# ╔═╡ 4cb153af-c1a0-47f5-a3d1-dfaa72f1c8b5
cell_type_legend # rename variable

# ╔═╡ 4c1a3c3b-f716-4437-8da7-0539673755e5
DATA_PATH = joinpath(dirname(@__DIR__),"data")

# ╔═╡ d40772a9-51a3-494a-b7f2-0bedcd2a2dc3
begin
	panel_df = CSV.read(joinpath(DATA_PATH,"panel.csv"), DataFrame);
	JLD2.@load joinpath(@__DIR__,"data","targets.jld2") targets
end

# ╔═╡ 43d75174-d346-40cf-ab10-40c50adf5b4f
begin
update_save
@bind hsv PlutoUI.combine() do Child
	md"""
	H $(Child(Slider(0:30:360))
	) S $(Child(Slider(range(1,0,21)))
	) V $(Child(Slider(range(1,0,21)))
	) 
	
	add target $(Child(Select(vcat([""],sort(targets))))
	) 
	"""
end
end

# ╔═╡ 85c3cc24-b7c0-41cf-a32a-9094d00d2304
IMG_PATH = joinpath(DATA_PATH,"img")

# ╔═╡ 812948b4-0efe-499c-be2f-3a1a26002ff0
FIG_PATH = joinpath(@__DIR__,"plots","images")

# ╔═╡ d6fbd6f8-fee6-4abd-b8e6-912e84f39954
selected_image = splitext(patient_image_dict[selected_patient])[1]

# ╔═╡ 434a387a-54ab-4147-89d2-e7fed1d7b0f4
	cell_coloring = Dict{Int,RGB{Float64}}(
		c=> if ci in highlighted
			highlight
		else
			get(cols, get(labels, ci, 0) ,unknown_color)  
		end
			for (c,ci) in eachrow(df_clustered[df_clustered.Image .== selected_image*".tiff", [:Object,:cluster_index]]))

# ╔═╡ bf024f5d-2ff4-42d4-a2fd-99f4f83df968
md"""
## Ref blocks for evaluation control
"""

# ╔═╡ e4bb254e-86ef-4c25-9219-66acb278efc9
begin # initialize csdict_saved
if in(readdir(joinpath(@__DIR__,"data")))("colorschemes.jld2") # if the file exists
 	JLD2.@load joinpath(@__DIR__,"data","colorschemes.jld2") csdict_saved
	else # initialize the file
	csdict_saved = Dict{String, Dict{String, HSV{Float32}}}()
	JLD2.@save joinpath(@__DIR__,"data","colorschemes.jld2") csdict_saved
end
available_keys = Ref{Vector{String}}([])
available_keys[] = collect(keys(csdict_saved));
end

# ╔═╡ a629f092-4613-4a67-a366-910b7184302c
@bind load_delete PlutoUI.combine() do Child
	md"""
	load: $(Child(Select(sort(available_keys[]) ))))
	delete: $(Child(Select(sort(available_keys[] ))))
	"""
end

# ╔═╡ e90a2701-a5e5-4eea-9f44-28c527a348d3
md"""
### saving
"""

# ╔═╡ c9e58092-bc6f-4c09-9f8c-5b3ba0339b62
load_delete_ref = ["",""]

# ╔═╡ 437dd0a9-3d9b-49d1-ac52-9a205bca38fa
begin
	load_delete_ref[1] = load_delete[1]
	load_delete_ref[2] = load_delete[2]
end

# ╔═╡ 5e261569-b068-4f8e-aeab-809f389759d5
color_collection = Dict{String,HSV{Float32}}() # rerun to clear

# ╔═╡ 1ada4235-c502-46d6-8861-8ae7b75fff68
begin 
	update_save # update the JLD2 with the current color scheme@
	default_name = load_delete_ref[1]
	if load_delete_ref[1] != "" 
		# first delete all elements of color_collection
		for k in keys(color_collection)
			delete!(color_collection,k)
		end
		# add all new keys to color collection
		merge!(color_collection, csdict_saved[load_delete_ref[1]])
	end
	load_delete_ref[2] != "" && delete!(csdict_saved, load_delete_ref[2]) # cs_delete is to to have a ref boundary
	JLD2.@save joinpath(@__DIR__,"data","colorschemes.jld2") csdict_saved
	available_keys[] = collect(keys(csdict_saved));
	color_update2 = true
end

# ╔═╡ 31c82837-417a-4da4-b713-448baaf377ea
@bind colorscheme_name PlutoUI.combine() do Child
	md"""
	colorscheme name $(Child(TextField(default=default_name)))
	"""
end

# ╔═╡ 8f3f5f9f-f815-4fef-9027-6de431ec9f2d
color_scheme = ["",HSV(1,1,0)]; # initialize a mutable container for controlling evaluation.

# ╔═╡ 2199af3d-9756-428b-8ecd-18fb34013abb
begin 
	color_update
	color_update2
	first(color_scheme) != "" && push!(color_collection, Pair(color_scheme...)) 
	included_markers = collect(keys(color_collection))
	cc = color_collection
end

# ╔═╡ f7858065-6c79-4636-af1d-5664d73cab85
begin 
update_save # set to default if we have loaded colors
@bind rmmarker PlutoUI.combine() do Child
	md"""
	drop target $(Child(Select(vcat([" "], sort(included_markers))))
	)
	"""
end
end

#@bind rmmarker Select(vcat([" "], sort(included_markers)))

# ╔═╡ 6fba6f0b-cce8-4a2a-a970-d6a5c9c2c2f9
begin
	H, S, V, marker = hsv
	color_scheme[2] = HSV(H,S,V)
	color_scheme[1] = marker
	delete!(color_collection, rmmarker[1])
	Pair(color_scheme...)
end

# ╔═╡ c375d2c7-a5c5-43d2-9de5-44754c17e1ac
function celltype_count_dict(cluster_vec, labels)
	# had trouble maintaining the order of the cells
	# returns a dict celltype => counts
	type_names = unique(collect(values(labels)))
	out = Dict(t => 0 for t in type_names)
	for c in cluster_vec
		out[labels[c]] += 1
	end
	return out
end

# ╔═╡ 9ad2fe2b-45bc-46bc-869f-51b1c051d705
cs = ["",Dict(""=>HSV(0,1,1))] # (name, Dict)

# ╔═╡ 3b08116a-e5f8-4e6b-a077-1d15c93bfa5e
begin
save_button	
csdict_saved[cs[1]] = deepcopy(cs[2]) # cs_add is to to have a ref boundary
JLD2.@save joinpath(@__DIR__,"data","colorschemes.jld2") csdict_saved
available_keys[] = collect(keys(csdict_saved));
end

# ╔═╡ bd0f192c-f2e3-4d55-9e15-2d1ce6f37030
let 
	cs[1] = colorscheme_name[1]
	cs[2] = color_collection
end

# ╔═╡ d2478db9-37a4-40cd-ab1d-7bef00352740
"""
generate the normalized bar plots for the cell type distribution across the images.
"""
function countplot_plot( highlighted_counts, total_counts,# list of celltypes
		image_order; # list of images
	image_names = img -> img, #get(image_patient_key,img,"(..."*img[end-17:end-5]*")"),
	highlight)
	colorfun = [RGB(.8), highlight] 
	
	hcounts = map(x->highlighted_counts[x], image_order)
	tcounts = map(x->total_counts[x], image_order) 
		
	# generate a matrix of cumulative probability
	cumprob_mat = permutedims(hcat(tcounts,hcounts),(2,1))
	
	fig2 = CairoMakie.Figure(size = (1400,600))
	ax2 = CairoMakie.Axis(fig2[1,1]; 
		aspect = 3, 
		xticks = (axes(cumprob_mat,2), map(image_names,image_order)), 
		xticklabelrotation = pi/2,
		xticklabelsize = 10);
	CairoMakie.xlims!(ax2, [.5,length(image_order)+.5])
	#CairoMakie.ylims!(ax2,0,1)

	for jj in axes(cumprob_mat,1)
		CairoMakie.barplot!(ax2,axes(cumprob_mat,2),cumprob_mat[jj,:], 
			color = colorfun[jj], gap = 0.0)
	end
	CairoMakie.vlines!(ax2,collect(axes(cumprob_mat,2)) .+ .5, color=RGB(1), linewidth = .5)
	
	# CairoMakie.save(joinpath(@__DIR__,"plots","cell_numbers.pdf"),fig2)
	fig2
end

# ╔═╡ 3eaceb5b-1d1d-4302-b542-118f79092527
 countplot_plot( highlighted_counts, total_counts, # list of celltypes
		collect(keys(highlighted_counts))[sortperm([
			-highlighted_counts[x] for x in collect(keys(highlighted_counts))])];highlight)

# ╔═╡ 12e96ec3-00cf-4f83-92b7-0ae10256d711
"""
generate the normalized bar plots for the cell type distribution across the images.
"""
function prob_plot(image_count_dict, 
		celltype_order, # list of celltypes
		image_order; # list of images
	image_names = img -> img, #get(image_patient_key,img,"(..."*img[end-17:end-5]*")"),
	colorfun = cell -> cols[cell])

		
	# generate a matrix of cumulative probability
	cumprob_mat = stack(
		cumsum(LinearAlgebra.normalize([image_count_dict[img][typ] for typ in celltype_order],1)) 
		for img in image_order) 
	
	fig2 = CairoMakie.Figure(size = (1400,600))
	ax2 = CairoMakie.Axis(fig2[1,1]; 
		aspect = 3, 
		xticks = (axes(cumprob_mat,2), map(image_names,image_order)), 
		xticklabelrotation = pi/2,
		xlabel = "Patient",
		ylabel = "Celltype fraction",
		xticklabelsize = 12);
	CairoMakie.xlims!(ax2, [.5,length(image_order)+.5])
	CairoMakie.ylims!(ax2,0,1)

	for jj in reverse(axes(cumprob_mat,1))
		CairoMakie.barplot!(ax2,axes(cumprob_mat,2),cumprob_mat[jj,:], 
			color = colorfun(celltype_order[jj]), gap = 0.0)
	end
	CairoMakie.vlines!(ax2,collect(axes(cumprob_mat,2)) .+ .5, color=RGB(1), linewidth = .5)
	
	fig2
end

# ╔═╡ 086e88f3-6505-4ef9-9b1e-1285290fa515
"""
helper function to convert the keys into an indexing scheme for fast colorizing
"""
function convert_to_sparse(color_scheme, panel_df)
	targets = keys(color_scheme)
	idx = Int[findfirst( ==(t), panel_df.name) for t in targets]
	colors = RGB{Float64}[RGB{Float64}(color_scheme[t]) for t in targets]
	return (idx,colors)
 end

# ╔═╡ 3a045540-2f87-49ec-9e46-5d446f61cc1d
"""
complicated function with the sole purporse of drawing the color wheel using luxor
"""
function draw_scale(hue,saturation; nh = 360, ns = 100)
    maxr = saturation + (1 / ns / 2 + .01)
    minr = saturation - (1 / ns / 2)
    maxθ = 2*pi*(hue/360 + (1 / nh / 2) + .001)
    minθ = 2*pi*(hue/360 - (1 / nh / 2))
    hsv = RGB(HSV(360-hue,saturation,1))
    sethue(red(hsv),green(hsv),blue(hsv))
    setline(0.1)
    arc(Point(0,0), maxr,minθ,maxθ)
    carc(Point(0,0), minr,maxθ,minθ)
    fillstroke()
end

# ╔═╡ 34f5a972-1aba-4dc3-8fa2-4de502f5f26e
@svg begin
    setline(10)
    sethue("purple")
    scale(.9)
    gsave()
    scale(100)
    for ii in 1:60
        for jj in range(0,1,20)
            draw_scale(ii*6,jj;nh = 60,ns = 20)
        end
    end
    grestore()
	gsave()
    sethue("black")
    setline(1.0)
    #for ii in 0:10:350
    #    θ = -(ii/360)*2*pi
    #    line(Point(100*cos(θ), 100*sin(θ)),Point(105*cos(θ), 105*sin(θ)); action = :stroke)
    #    text(string(ii), Point(110*cos(θ), 110*sin(θ)), angle=θ, valign = :middle)
    #end
	grestore()
	sethue("black")
	setline(3)
	circle(Point(S*100*cospi(H/180), -S*100*sinpi(H/180)),4, action = :stroke)
end 200 200 FIG_PATH*"/.tmp_colorwheel"

# ╔═╡ cdf7bade-38dc-45da-b685-628009c8e945
"""
Generates the color images using the loaded images and color collection
"""
function colorize(tiff_array, color_scheme; 
	panel_df = panel_df, value_scale = 1)
	idx, clrs = convert_to_sparse(color_scheme, panel_df)
	img = fill(RGB(0), size(tiff_array,2), size(tiff_array,3)) # allocate an array of RGB
	for I in CartesianIndices(img)
		pixel = RGB(0)
		for (i,c) in zip(idx, clrs)
			pixel += tiff_array[i,I] * c * value_scale # for each color c, add it to the pixel
		end
		img[I] = RGB{Float64}(clamp(pixel.r,0.0,1.0), clamp(pixel.g,0.0,1.0), clamp(pixel.b,0.0,1.0)) # clamp dynamic range
	end
	return img
end

# ╔═╡ 9f9be021-8ad3-409a-bdb8-eac51b64a504
getimage_tiff(name) = channelview(TiffImages.load(joinpath(DATA_PATH,"img",name*".tiff"), verbose = false))

# ╔═╡ 56ad0689-2552-4098-961c-1a64ad326783
@memoize function image_channel_scales(name; reference_quantile = 0.95)
	nan_remove(x) = isinf(x) ? (reference_quantile*255) : 
		clamp(Float64(x), 0.5, 255)
	data = getimage_tiff(name)
	return [min(
		(reference_quantile) /
			nan_remove(quantile(vec(view(data, :, :, k)), reference_quantile)),1.0) for k in axes(data,3)]
end

# ╔═╡ 5891fc07-1741-4268-9ee5-d0f803ce9e7b
"""
input 'data' is an array [x,y,channel]
returns quantile scaled between 0 and 1 with [channel, x,y]
reference_quantile is a parameter for the normalization method
"""
function quantile_scale(data; reference_quantile = 0.95, name = nothing)
	out = similar(data)
	nan_remove(x) = isinf(x) ? (reference_quantile*255) : 
		clamp(Float64(x), 0.5, 255)
	vec_scales::Vector{Float64} = image_channel_scales(name; reference_quantile = 0.95)
	for k in axes(data,3)
		#scale_factor = (reference_quantile) /
		#	nan_remove(quantile(vec(view(data, :, :, k)), reference_quantile))
		scale_factor = vec_scales[k]  # min(scale_factor,1)
		for j in axes(data,2)
			for i in axes(data,1)		
				out[i,j,k] = clamp( scale_factor * data[i, j, k], 0, 1)
			end
		end
	end
	return permutedims(out,(3,1,2))
end

# ╔═╡ 89548f46-bb68-4ee5-bedb-b310c916f359
"""
loads an image by name and returns the colorized version
"""
function color_img(name, color_scheme; kwds...)
	arr = quantile_scale(getimage_tiff(name); name)
	colorize(arr, color_scheme; kwds...)
end

# ╔═╡ c5f01676-388f-4bca-bf85-af3ced1bcc2d
getmask_tiff(name) = Matrix{Int}(rawview(channelview(TiffImages.load(
        joinpath(DATA_PATH,"masks",name*".tiff"), verbose = false))))

# ╔═╡ 9dd8c509-d818-45d8-84d9-52e5b95f0dcd
"""
Save images in the folder 
"""
function color_img_save(name, img, color_scheme; value_scale = 1.0, window = nothing)
	picture = color_img(img, color_scheme; value_scale)
    s = size(picture)
    if isnothing(window)
        window = (1:s[1],1:s[2])
    end
	Drawing(length(window[2]),length(window[1]), 
        # where we save the image
        joinpath(FIG_PATH, img*name))
        # begin image
    gsave()
    Luxor.transform([0 1 1 0 0 0])
    gsave()
        Luxor.transform([0 1 1 0 0 0])
        placeimage(picture[window[1], window[2]], Luxor.Point(0, 0))
        grestore()
	finish()
end

# ╔═╡ 6858f726-fde2-4273-a1c5-5268617a53d6
"""
draw legends and return the bounding box coordinates
"""
function draw_legends(pairs...)
    vec_wrap(x::Vector) = x
    vec_wrap(x::AbstractString) = [x]
    loremipsum = [vec_wrap(jj) for (jj,_) in pairs]
    colors = [ii for (_,ii) in pairs]
    fontface("Helvetica")
    fontsize(20)
    # closure that captures the widest line
    _counter() = (w = 0; (n) -> w = max(w,n))
    counter = _counter()

    setopacity(0)
    h = textbox(map(filter(!isempty,loremipsum)) do x
        join(x, " + ")
    end,
        O + (10, 0),
        leading = 30,
        linefunc = (lnumber, str, pt, h) -> begin
            sethue(colors[lnumber])
            counter(textextents(str)[3])
        end)
    setopacity(1.0)
    sethue(RGBA(0,0,0))
	background(RGBA(0,0,0))
	bb = BoundingBox(Luxor.box(O, Luxor.Point(O.x + counter(0) + 20, h.y - 15); vertices = true))
    #box(bb,10 ;action=:fill)
    setopacity(1.0)
    textbox(map(filter(!isempty,loremipsum)) do x
        join(x, " + ")
    end,
        O + (10, 0),
        leading = 30,
        linefunc = (lnumber, str, pt, h) -> begin
        sethue(colors[lnumber])
        end)
	return bb
end

# ╔═╡ 5b15bb55-6b81-47e3-beef-a292bd20562c
"""
Save the legend
"""
function color_scheme_legend(name, color_scheme; value_scale = 1.0)
	if isempty(color_scheme)
		color_scheme = Dict("no colors"=>RGB(1,1,1))
	end
	Drawing(100, 100, "/tmp/emptydrawing.png")
	bb = draw_legends(pairs(color_scheme)...)
	finish()
	Luxor.@svg begin
		translate(- getfield(bb.corner2,:x)/2, - getfield(bb.corner2,:y)/2)
		draw_legends(pairs(color_scheme)...)
	end getfield(bb.corner2,:x) getfield(bb.corner2,:y) "/tmp/tmp_legend.svg"
end

# ╔═╡ 4e6989c4-72b6-4b54-ae6e-0a035e6cc8e1
marker_legend = color_scheme_legend("scheme", cc; value_scale = 1.0) 

# ╔═╡ ed85f170-0bc4-4029-9652-b20f27000894
marker_legend

# ╔═╡ 488fb0d5-f61a-43c5-a85b-8a7a3542e8ae
roi_figure = 
	let 
	if "cells" in cellbits 
		out = map(i -> ifelse(i == 0, RGB(0.0), get(cell_coloring,i, unknown_color)), getmask_tiff(selected_image))
		name = "cells"
	else	
		out = color_img(selected_image, cc; value_scale = value_scale)
		name = cs[1]
	end
		
	if "save image" in cellbits  
		save(joinpath(@__DIR__,"plots","images",selected_patient * "_" * name * ".png"), out)
		if "cells" in cellbits 
			save(joinpath(@__DIR__,"plots","images","0_cell_type_legend.pdf"), cell_type_legend)
		else	
			save(joinpath(@__DIR__,"plots","images","00_legend_$(cs[1]).svg"), marker_legend)
		end
	end
		
	out
end;

# ╔═╡ 06f9feb8-6ef1-4c04-af42-c75344f6509e
@bind box RectangleDrawing(roi_figure)

# ╔═╡ d0e5885f-44ed-4223-8fcd-1694ad13b09e
try FullImage(repeat(map(i -> ifelse(i == 0, RGB(0.0), get(cell_coloring,i, unknown_color)), getmask_tiff(selected_image)[get_bounds(box...)...]),inner = (7,7)))
catch
	md"""
	choose region
	"""
end	

# ╔═╡ c430663a-d72e-426e-9f66-17dbb6964f19
cc

# ╔═╡ 6ca8e213-3b45-44fc-989a-6a938b1ff785
md"""
# Drawing Single Cell boundaries 
"""

# ╔═╡ c28e4ba4-4671-4cc2-a3eb-5ee8b8d0aee6

function addpoints!(bd, v1, i1, j1, v2, i2, j2)
    imid = (i1+i2) / 2
    jmid = (j1+j2) / 2
    idif = (i1-i2) / 2
    jdif = (j1-j2) / 2

    push!(bd[v1], 
        (Luxor.Point(imid - jdif, jmid + idif), Luxor.Point(imid + jdif, jmid - idif)) )
    push!(bd[v2], 
        (Luxor.Point(imid + jdif, jmid - idif), Luxor.Point(imid - jdif, jmid + idif)) )
end


# ╔═╡ 69fec845-702e-409d-9b61-c118706d26a5
"""
this function scans through a mask tiff, that is an array of integers associated with the cell object and generates a dictionary that contains the points neccesary to draw the polygon for the cell boundary
"""
function cellboundaries(mask)
    (width,height) = size(mask)
    # returns an object
    bndry_dict = Dict(ii => Vector{Tuple{Luxor.Point,Luxor.Point}}() for ii in unique(mask))
    #draw verticle lines
    for ii in 1:width
        val = zero(mask[ii,1])
        for jj in 1:height
            new = mask[ii,jj]
            if new != val
                addpoints!(bndry_dict, val, ii, jj-1, new, ii, jj)
                val = new
            end
        end
        addpoints!(bndry_dict, val, ii, height, zero(val), ii, height + 1)
    end
    for jj in 1:height
        val = zero(mask[1,jj])
        for ii in 1:width
            new = mask[ii,jj]
            if new != val
                addpoints!(bndry_dict, val, ii-1, jj, new, ii, jj)
                val = new 
            end
        end
        addpoints!(bndry_dict, val, width, jj, zero(val), width+1, jj)
    end
    return bndry_dict
end

# ╔═╡ 1569e4dd-ee81-4824-8613-96126f51341b
function movecloser!(vec,r)
    for (ii, v0,v1,v2) in Iterators.zip(eachindex(vec), circshift(vec,1), vec, circshift(vec,-1))
        dif = v2 - v0
        rot_dif = Luxor.Point(-sign(dif.y), sign(dif.x))
        vec[ii] = v1 + r * rot_dif
    end
end


# ╔═╡ 5068f069-1811-4541-8d3e-983cf502c136
#rads(pt) = atan(pt.y, pt.x)
crossz(pt1,pt2) = pt1.y*pt2.x - pt1.x * pt2.y

# ╔═╡ d9d835d7-7b81-473a-975c-767f0708d184

"""
this function actually draws the cell boundaries in a particular hue within a luxor image.
There is an entry associated with the region of pixels that have not been segemented

"""
function draw_singlecell_boundaries(vd)
    ii = 1
    outs = Vector{Vector{Luxor.Point}}()
    out = Vector{Luxor.Point}()
    while true
        (v1,v2) = popat!(vd, ii)
        push!(out,v1)
        isempty(vd) ? break : nothing
        ii_list = findall(x-> x[1] == v2, vd)
        if isempty(ii_list) # start over if you can't find your way
            push!(outs,out)
            out = Vector{Luxor.Point}()
            ii = 1
        else
            ii = argmin( ii -> crossz(v1-v2, vd[ii][2] - v2), ii_list)
        end
    end
    push!(outs,out)
    sethue(RGB(.8))
    for out in outs
        movecloser!(out,.2)
        Luxor.poly(out; action = :stroke, close = true)
    end
end

# ╔═╡ 1957d5e3-c9a1-4229-9f42-b502c21ddfb0
"""
draws all the boundaries in the dictionary
"""
function draw_boundarydict(bd)
    for (val, vec) in pairs(bd)
        if val > 0
            draw_singlecell_boundaries(vec)
        end
    end
end


# ╔═╡ 1aa5bd45-efb8-4fee-ba87-874cf3b4b2b9
md"""
# Drawing full color figures with masks
"""

# ╔═╡ f3123e19-b4c9-4ed7-b59c-032e49c420ba
FIG_PATH

# ╔═╡ e3d16459-06bd-4812-8e2f-e2830707d731
typeof(range(1,2))

# ╔═╡ 3cc89a3e-5405-496e-94b4-9a5165558915
function color_figure(img; mask_dict = nothing, 
    window = nothing, # optional parameter for limiting the size of the image
	name = "untitled_$(rand(0:9,4)...)",
    figpath = FIG_PATH, colorscale = 0.7, scalefactor = 5)
	
    s = size(img)
    if isnothing(window)
        window = (1:s[1],1:s[2])
	else
		window = (max(1,minimum(window[1])):min(s[1],maximum(window[1])),max(1,minimum(window[2])):min(s[2],maximum(window[2])))
    end

    a = Drawing(scalefactor*length(window[2]),scalefactor*length(window[1])..., 
        # where we save the image
        joinpath(figpath, name*".svg"))
        # begin image
	scale(scalefactor)
    gsave()
    Luxor.transform([0 1 1 0 0 0])
        gsave()
            Luxor.transform([0 1 1 0 0 0])
            placeimage(img[window[1], window[2]], Luxor.Point(0, 0))
        grestore()
        gsave()
        if !isnothing(mask_dict)
            setline(1.0)
            translate(Luxor.Point(0.5 - minimum(window[1]), 0.5 - minimum(window[2])))
            draw_boundarydict(mask_dict)
        end
        grestore()
    grestore()
    finish()
	return a
end

# ╔═╡ 02e89f20-f528-4abe-8596-f3aced0bdd48
try color_figure(color_img(selected_image, cc; value_scale = value_scale); 
	mask_dict = 
		ifelse("masks" in maskbits, cellboundaries(getmask_tiff(selected_image)),nothing),
    window = get_bounds(box...), # optional parameter for limiting the size of the image
	name = ".tmp_widow",
    figpath = FIG_PATH, colorscale = value_scale, scalefactor = 7)
catch
	md"""
	choose region
	"""
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Base64 = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
ImageCore = "a09fc81d-aa75-5fe9-8630-4744c3626534"
ImageIO = "82e4d734-157c-48bb-816b-45c225c6df19"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
Memoization = "6fafb56a-5788-4b4e-91ca-c0cea6611c73"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
TiffImages = "731e570b-9d59-4bfa-96dc-6df516fadf69"

[compat]
CSV = "~0.10.14"
CairoMakie = "~0.12.13"
DataFrames = "~1.7.0"
FileIO = "~1.16.3"
HypertextLiteral = "~0.9.5"
ImageCore = "~0.10.2"
ImageIO = "~0.6.8"
JLD2 = "~0.5.5"
Luxor = "~4.1.0"
Memoization = "~0.2.1"
PlutoUI = "~0.7.60"
TiffImages = "~0.10.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "094ca28ff59a0c8fc5d7ce72e3818f20436a703c"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "014bc22d6c400a7703c0f5dc1fdc302440cf88be"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.4"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "6c834533dc1fabd820c1db03c839bf97e45a3fab"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.14"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "7b6ad8c35f4bc3bca8eb78127c8b99719506a5fb"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.0"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "2b04b60ed9d3e977f93e34952971b608c34b3401"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.12.13"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "PrecompileTools", "Random"]
git-tree-sha1 = "668bb97ea6df5e654e6288d87d2243591fe68665"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "d7477ecdafb813ddee2ae727afa94e9dcb5f3fb0"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.112"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.Extents]]
git-tree-sha1 = "81023caa0021a41712685887db1fc03db26f41f5"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.4"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "7878ff7172a8e6beedd1dea14bd27c3c6340d361"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.22"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "2493cdfd0740015955a8e46de4ef28f49460d8bc"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.3"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "59107c179a586f0fe667024c5eb7033e81333271"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.2"

[[deps.GeoInterface]]
deps = ["Extents", "GeoFormatTypes"]
git-tree-sha1 = "2f6fce56cdb8373637a6614e14a5768a88450de2"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.7"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "fc713f007cff99ff9e50accba6373624ddd33588"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "b2a7eaa169c13f5bcae8131a83bc30eff8f71be0"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.2"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "8e125d40cae3a9f4276cdfeb4fcdb1828888a4b3"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.17"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

    [deps.IntervalArithmetic.weakdeps]
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "Requires", "TranscodingStreams"]
git-tree-sha1 = "aeab5c68eb2cf326619bf71235d8f4561c62fe22"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.5.5"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "f389674c99bfcde17dc57454011aa44d5a260a40"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "PolygonAlgorithms", "PrecompileTools", "Random", "Rsvg"]
git-tree-sha1 = "134570038473304d709de27384621bd0810d23fa"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "4.1.0"
weakdeps = ["LaTeXStrings", "MathTeXEngine"]

    [deps.Luxor.extensions]
    LuxorExtLatex = ["LaTeXStrings", "MathTeXEngine"]

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "50ebda951efaa11b6db0413c1128726b8eab3bf0"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.21.13"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "4604f03e5b057e8e62a95a44929cafc9585b0fe9"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.8.9"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "e1641f32ae592e415e3dbae7f4a188b5316d4b62"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.1"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Memoization]]
deps = ["MacroTools"]
git-tree-sha1 = "073f080e733bc6697411901224ed4fd15fefaffa"
uuid = "6fafb56a-5788-4b4e-91ca-c0cea6611c73"
version = "0.2.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PolygonAlgorithms]]
git-tree-sha1 = "a5ded6396172cff3bacdd1354d190b93cb667c4b"
uuid = "32a0d02f-32d9-4438-b5ed-3a2932b48f96"
version = "0.2.0"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "98ca7c29edd6fc79cd74c61accb7010a4e7aee33"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.6.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ff11acffdb082493657550959d4feb4b6149e73a"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "bc7fd5c91041f44636b2c134041f7e5263ce58ae"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "86e7731be08b12fa5e741f719603ae740e16b666"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.10+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "7dfa0fd9c783d3d0cc43ea1af53d69ba45c447df"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═a112be8c-2245-11ef-2c0f-1d00ed8b2f19
# ╠═cbe74fb5-972f-4127-a6c7-0c50cf8522cb
# ╠═28214a08-028a-415c-b294-8a290931e8eb
# ╠═9a87684f-a7e6-4301-af1c-06166546f344
# ╠═9fcdd12a-cd55-47b4-a98a-6d1c84f4b436
# ╟─7bce7789-58d5-49e7-b0cc-f38f2829f0db
# ╠═73279267-ba04-47da-9bba-81a1bc26269c
# ╟─34f5a972-1aba-4dc3-8fa2-4de502f5f26e
# ╠═43d75174-d346-40cf-ab10-40c50adf5b4f
# ╟─6fba6f0b-cce8-4a2a-a970-d6a5c9c2c2f9
# ╟─f7858065-6c79-4636-af1d-5664d73cab85
# ╟─9e355daf-ada9-4519-b186-2512a63a6747
# ╟─2199af3d-9756-428b-8ecd-18fb34013abb
# ╟─714ed9c9-7307-482c-b9dc-5c5027aae3e0
# ╟─a629f092-4613-4a67-a366-910b7184302c
# ╟─f501dad3-6a8f-43a6-9931-75cdc7451769
# ╟─ffd766c6-e030-484b-a769-0bd3bad19fa4
# ╟─31c82837-417a-4da4-b713-448baaf377ea
# ╟─50892ac6-3480-4f1b-9f9b-bf60753989c4
# ╟─a79670ce-9a64-41a8-939f-4cb4e6a84dfe
# ╠═39f58298-320a-4e7e-a170-7516d85149d4
# ╟─4e6989c4-72b6-4b54-ae6e-0a035e6cc8e1
# ╟─06f9feb8-6ef1-4c04-af42-c75344f6509e
# ╟─c0170387-a1d1-4f17-a41c-ae3875c2034a
# ╟─8272ba6b-b303-43ad-8b3f-a96e28e387e0
# ╠═0aecf292-b442-4b9a-b690-06e061a6af50
# ╟─ba317e1a-c67c-4f5a-b62a-56d9c0f0df86
# ╟─ce169cd6-91ab-4efb-9915-88d09f36038c
# ╠═ed85f170-0bc4-4029-9652-b20f27000894
# ╟─02e89f20-f528-4abe-8596-f3aced0bdd48
# ╟─d0e5885f-44ed-4223-8fcd-1694ad13b09e
# ╟─4cb153af-c1a0-47f5-a3d1-dfaa72f1c8b5
# ╟─0c29afcf-c687-4776-b870-84583727592c
# ╟─da0cc366-5471-4f93-bc6d-74739910cffb
# ╠═3eaceb5b-1d1d-4302-b542-118f79092527
# ╠═7d73f953-4365-433d-b7f4-0716c79805f0
# ╠═b558f09a-9468-4b49-915b-a614ee009618
# ╠═f799657d-681e-42a2-9261-f1ab1c11421f
# ╟─4d9deb1e-77e6-44ca-9d7d-fbb7101e4f39
# ╠═d40772a9-51a3-494a-b7f2-0bedcd2a2dc3
# ╠═718eb743-7b7c-46dd-890b-18a51240db2e
# ╠═86841954-e590-45a9-936b-7bd653109f4c
# ╠═08358c11-cf3f-49a1-bde9-f143899c1557
# ╠═3daf7f1c-4687-468d-80f8-03ede23b39e8
# ╟─952f5688-ba41-49ec-b5d7-f3e04cf9b4bc
# ╟─6eac3e81-d949-4a59-ae77-0c16a7f29aaa
# ╠═339eb3c7-40de-4a5d-82b2-c571cf421288
# ╠═434a387a-54ab-4147-89d2-e7fed1d7b0f4
# ╠═85c3cc24-b7c0-41cf-a32a-9094d00d2304
# ╠═4c1a3c3b-f716-4437-8da7-0539673755e5
# ╠═812948b4-0efe-499c-be2f-3a1a26002ff0
# ╠═d6fbd6f8-fee6-4abd-b8e6-912e84f39954
# ╠═488fb0d5-f61a-43c5-a85b-8a7a3542e8ae
# ╟─bf024f5d-2ff4-42d4-a2fd-99f4f83df968
# ╠═e4bb254e-86ef-4c25-9219-66acb278efc9
# ╠═3b08116a-e5f8-4e6b-a077-1d15c93bfa5e
# ╠═1ada4235-c502-46d6-8861-8ae7b75fff68
# ╟─e90a2701-a5e5-4eea-9f44-28c527a348d3
# ╠═c9e58092-bc6f-4c09-9f8c-5b3ba0339b62
# ╠═437dd0a9-3d9b-49d1-ac52-9a205bca38fa
# ╠═5e261569-b068-4f8e-aeab-809f389759d5
# ╠═8f3f5f9f-f815-4fef-9027-6de431ec9f2d
# ╠═c375d2c7-a5c5-43d2-9de5-44754c17e1ac
# ╠═bd0f192c-f2e3-4d55-9e15-2d1ce6f37030
# ╠═9ad2fe2b-45bc-46bc-869f-51b1c051d705
# ╠═d2478db9-37a4-40cd-ab1d-7bef00352740
# ╠═12e96ec3-00cf-4f83-92b7-0ae10256d711
# ╠═086e88f3-6505-4ef9-9b1e-1285290fa515
# ╠═3a045540-2f87-49ec-9e46-5d446f61cc1d
# ╠═cdf7bade-38dc-45da-b685-628009c8e945
# ╠═5891fc07-1741-4268-9ee5-d0f803ce9e7b
# ╠═56ad0689-2552-4098-961c-1a64ad326783
# ╠═89548f46-bb68-4ee5-bedb-b310c916f359
# ╠═9f9be021-8ad3-409a-bdb8-eac51b64a504
# ╠═c5f01676-388f-4bca-bf85-af3ced1bcc2d
# ╠═9dd8c509-d818-45d8-84d9-52e5b95f0dcd
# ╠═6858f726-fde2-4273-a1c5-5268617a53d6
# ╠═5b15bb55-6b81-47e3-beef-a292bd20562c
# ╠═c430663a-d72e-426e-9f66-17dbb6964f19
# ╟─6ca8e213-3b45-44fc-989a-6a938b1ff785
# ╠═c28e4ba4-4671-4cc2-a3eb-5ee8b8d0aee6
# ╠═69fec845-702e-409d-9b61-c118706d26a5
# ╠═d9d835d7-7b81-473a-975c-767f0708d184
# ╟─1569e4dd-ee81-4824-8613-96126f51341b
# ╟─1957d5e3-c9a1-4229-9f42-b502c21ddfb0
# ╟─5068f069-1811-4541-8d3e-983cf502c136
# ╟─1aa5bd45-efb8-4fee-ba87-874cf3b4b2b9
# ╠═f3123e19-b4c9-4ed7-b59c-032e49c420ba
# ╠═e3d16459-06bd-4812-8e2f-e2830707d731
# ╠═3cc89a3e-5405-496e-94b4-9a5165558915
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
