#!/usr/bin/env python

from __future__ import print_function
import os, sys, argparse, json, matplotlib.pyplot, seaborn

seaborn.set()
seaborn.set_style("white")

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#*******************************************************************************

class Plotter:
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        pass

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def draw_pile(self, title, pile):
        if ("data" not in pile or "begin" not in pile or "end" not in pile or
            "chimeric" not in pile or "repetitive" not in pile or
            "median" not in pile):
            eprint("[raven::Plotter::plot_pile] error: incomplete pile!")
            return

        figure, ax = matplotlib.pyplot.subplots(1, 1, figsize=(7.5, 5))

        scpb = seaborn.color_palette("Blues")
        scpr = seaborn.color_palette("Reds")

        x = range(len(pile["data"]))
        seaborn.despine()

        ax.plot(x, pile["data"], label="data", color=scpb[2])
        for slope in pile["chimeric"]:
            ax.axvline(slope, color=scpr[2], linestyle=":")
        for slope in pile["repetitive"]:
            ax.axvline(slope, color=scpr[3], linestyle=":")
        begin = pile["begin"]
        end = pile["end"]

        ax.axhline(int(pile["median"]), label="median", color=scpb[1], linestyle=":")
        ax.set_title(title)
        ax.set_ylim([0, int(pile["median"]) * 3])

        figure.text(0.5, 0.04, "base", ha="center")
        figure.text(0.04, 0.5, "coverage", va="center", rotation="vertical")
        matplotlib.pyplot.legend(loc="best")
        matplotlib.pyplot.savefig(str(title) + ".pdf", format='pdf', dpi=1200)
        matplotlib.pyplot.close(figure)

    def draw_assembly(self, title, component):
        if ("nodes" not in component or "edges" not in component):
            eprint("[raven::Plotter::plot_assembly] error: incomplete component!")
            return

        scpg = seaborn.cubehelix_palette(rot=-.4)
        scpr = seaborn.color_palette("Reds")

        matplotlib.pyplot.figure(figsize=(16,16), frameon=False)

        for edge in component["edges"]:
            x = component["nodes"][edge[0]]
            y = component["nodes"][edge[1]]
            if (x[2] == 1 or y[2] == 1):
                clr = scpr[3]
            else:
                clr = scpg[2] if edge[2] == 0 else scpg[1]
            matplotlib.pyplot.plot([x[0], y[0]], [x[1], y[1]], '-' if edge[2] == 0 else ':', c=clr)

        for node in component["nodes"]:
            x = component["nodes"][node]
            matplotlib.pyplot.plot(x[0], x[1], '.', c=scpg[4], markersize=(5 if x[3] == 1 else 15))

        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.yticks([])
        matplotlib.pyplot.axis('off')
        matplotlib.pyplot.savefig(title + '.pdf', format='pdf', dpi=1200)

    def run(self):
        try:
            d = open(self.path)
        except Exception:
            eprint("[raven::Plotter::run] error: unable to open file {}!".format(self.path))
            return

        try:
            data = json.load(d)
        except Exception:
            eprint("[raven::Plotter::run] error: file is not in JSON format!")
            return

        if (("piles" not in data or not data["piles"]) and
            ("assembly" not in data or not data["assembly"])):
            eprint("[raven::Plotter::run] error: incomplete input file!")
            return

        if ("piles" in data):
            for pile in data["piles"]:
                self.draw_pile(pile, data["piles"][pile])
        else:
            for component in data["assembly"]:
                self.draw_assembly(component, data["assembly"][component])

#*******************************************************************************

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Plot is a tool for drawing
        different stages of the raven assembler""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("path", help="""Input file in JSON format""")

    args = parser.parse_args()
    plotter = Plotter(args.path)

    with plotter:
        plotter.run()
