using Documenter

push!(LOAD_PATH, "../src/")

using ImageMetrics
import ImageMetrics

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    sitename = "Image Metrics",
    authors = "Éric Thiébaut and contributors",
    pages = ["index.md", "general.md", "bc2016.md"],
    format = Documenter.HTML(
        prettyurls = DEPLOYDOCS,
        mathengine = Documenter.KaTeX(
            Dict(
                :delimiters => [
                    Dict(:left => raw"$",   :right => raw"$",   display => false),
                    Dict(:left => raw"$$",  :right => raw"$$",  display => true),
                    Dict(:left => raw"\[",  :right => raw"\]",  display => true),
                ],
                :macros => Dict(raw"\Vh" => raw"\boldsymbol{h}",
                                raw"\Vt" => raw"\boldsymbol{t}",
                                raw"\Vx" => raw"\boldsymbol{x}",
                                raw"\Vy" => raw"\boldsymbol{y}",
                                raw"\Vz" => raw"\boldsymbol{z}",
                                raw"\Dist" => raw"\mathcal{D}",
                                raw"\Score" => raw"\mathcal{S}",
                                raw"\RR" => raw"\mathbb{R}",
                                raw"\MR" => raw"\mathbf{R}",
                                raw"\One" => raw"\boldsymbol{1}",
                                raw"\Zero" => raw"\boldsymbol{0}",
                                #raw"\Sign" => raw"\mathrm{sign}",
                                ),
            )
        ),
    ),
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/JMMC-OpenDev/ImageMetrics",
    )
end
