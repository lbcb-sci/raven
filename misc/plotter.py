#!/usr/bin/env python
import os, sys, argparse, json, seaborn
from matplotlib import pyplot

seaborn.set()
seaborn.set_style("white")
seaborn.despine()
scpb = seaborn.color_palette("Blues")
scpr = seaborn.color_palette("Reds")
scpg = seaborn.cubehelix_palette(rot=-.4)

class Plotter:
  def __init__(self, mode, path, type):
    self.mode = mode
    self.path = path
    self.type = type

  def DrawPile(self, pile):
    if ((self.type == "regular" and (pile["is_chimeric_"] or pile["is_repetitive_"])) or
        (self.type == "chimeric" and not pile["is_chimeric_"]) or
        (self.type == "repetitive" and not pile["is_repetitive_"])):
      return

    figure, ax = pyplot.subplots(1, 1, figsize = (7.5, 5))

    ax.plot(range(len(pile["data_"])), pile["data_"], label = "data", color = scpb[2])

    ax.axhline(int(pile["median_"]), label = "median", color = scpb[1], linestyle = ":")

    for slope in pile["chimeric_regions_"]:
      ax.axvline(slope["first"], color = scpr[2], linestyle = ":")
      ax.axvline(slope["second"], color = scpr[2], linestyle = ":")

    for slope in pile["repetitive_regions_"]:
      c = scpg[3] if slope["first"] & 1 else scpr[3]
      ax.axvline(slope["first"] >> 1, color = c, linestyle = ":")
      ax.axvline(slope["second"], color = c, linestyle = ":")

    ax.set_title(pile["id_"])
    ax.set_ylim([0, int(pile["median_"]) * 3])
    figure.text(0.5, 0.05, "base", ha = "center")
    figure.text(0.05, 0.5, "coverage", va = "center", rotation = "vertical")
    pyplot.legend(loc="best")
    pyplot.savefig(str(pile["id_"]) + ".png", format = 'png')
    pyplot.close(figure)

  def DrawGraph(self, title, graph):
    pyplot.figure(figsize = (16,16), frameon = False)

    for edge in graph["edges"]:
      x = graph["nodes"][edge[0]]
      y = graph["nodes"][edge[1]]
      if (x[2] == 1 or y[2] == 1):
        c = scpr[3]
      else:
        c = scpg[2] if edge[2] == 0 else scpg[1]
      pyplot.plot([x[0], y[0]], [x[1], y[1]], '-' if edge[2] == 0 else ':', color = c)

    for node in graph["nodes"]:
      x = graph["nodes"][node]
      pyplot.plot(x[0], x[1], '.', color = scpg[4], markersize = (5 if x[3] == 1 else 15))

    pyplot.xticks([])
    pyplot.yticks([])
    pyplot.axis('off')
    pyplot.savefig(title + '.pdf', format = 'pdf')

  def Run(self):
    try:
      f = open(self.path)
    except Exception:
      print("[raven::Plotter::Run] error: unable to open file {}".format(self.path))
      return
    try:
      data = json.load(f)
    except Exception:
      print("[raven::Plotter::Run] error: file is not in JSON format")
      return

    if (self.mode == "graph"):
      for component in data:
        self.DrawGraph(component, data[component])
    elif (self.mode == "pile"):
      for pile in data:
        self.DrawPile(data[pile])
    return

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
      description = "Plotter is a tool for drawing the assembly graph and pile-o-grams",
      formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("mode",
      help = "draw either the assembly [graph] or [pile]-o-grams")
  parser.add_argument("path",
      help = "input file in JSON format")
  parser.add_argument("-t", "--type", default = "all",
      help = "pile type selection [regular, chimeric, repetitive]")

  args = parser.parse_args()
  plotter = Plotter(args.mode, args.path, args.type)
  plotter.Run()
