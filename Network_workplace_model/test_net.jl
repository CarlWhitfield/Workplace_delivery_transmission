using Graphs
using LightGraphs
using GraphPlot
using Compose
using Cairo
using Colors

function make_and_plot_graph(title::String)



    g = Graphs.simple_graph(6,is_directed=false)
    lg = LightGraphs.SimpleGraph(6)
    ecount = 0
    while ecount < 10
        s = rand(1:6,1)[1]
        t = s
        while t == s
            t = rand(1:6,1)[1]
        end
        if !LightGraphs.has_edge(lg,s,t)
            LightGraphs.add_edge!(lg,s,t)
            Graphs.add_edge!(g,s,t)
            ecount += 1
        end
    end
    eind = edge_index.(Graphs.edges(g))
    colormap = [colorant"red",colorant"blue",colorant"yellow"]
    ecols = colormap[mod.(eind,3) .+ 1]

    draw(PNG(title,8cm,8cm),gplot(lg,edgestrokec=ecols,edgelabel=1:ne(lg),
                                  nodelabel=1:nv(lg)))
end
