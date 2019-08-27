#!/usr/bin/env python

from __future__ import print_function
import os, sys, argparse, json, matplotlib.pyplot, seaborn

seaborn.set()
seaborn.set_style("white")

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#*******************************************************************************

class Plotter:
    def __init__(self, data_path, ylimit, out_path):
        self.data_path = data_path
        self.ylimit = int(ylimit)
        self.out_path = out_path + "/"
        if (not os.path.isdir(self.out_path)):
            eprint("[rala::Plotter::__init__] error: invalid out directory!")
            sys.exit(1)

    def __enter__(self):
        pass

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    @staticmethod
    def plot_pile(pile, orientation, overlap_length, type, title, ylimit, ax):
        if ("data" not in pile or "begin" not in pile or "end" not in pile or
            "median" not in pile):
            eprint("[rala::Plotter::plot_pile] error: incomplete pile!")
            sys.exit(1);

        scpb = seaborn.color_palette("Blues")
        scpr = seaborn.color_palette("Reds")

        x = range(len(pile["data"]))
        ax.set_ylim([0, ylimit])
        seaborn.despine()

        if (orientation == 0):
            ax.plot(x, pile["data"], label="data", color=scpb[2])
            begin = pile["begin"]
            end = pile["end"]
        else:
            ax.plot(x, list(reversed(pile["data"])), label="data", color=scpb[2])
            begin = len(x) - pile["end"]
            end = len(x) - pile["begin"]

        ax.axhline(int(pile["median"]), label="median", color=scpb[1], linestyle=":")
        ax.set_title(title)

    def run(self):
        try:
            d = open(self.data_path)
        except Exception:
            eprint("[rala::Plotter::run] error: unable to open file {}!".format(self.data_path))
            sys.exit(1)

        try:
            data = json.load(d)
        except Exception:
            eprint("[rala::Plotter::run] error: file is not in JSON format!")
            sys.exit(1)

        if ("piles" not in data or not data["piles"]):
            eprint("[rala::Plotter::run] error: incomplete input file!")
            sys.exit(1)

        for pile in data["piles"]:
            figure, ax = matplotlib.pyplot.subplots(1, 1, figsize=(7.5, 5))

            Plotter.plot_pile(data["piles"][pile], 0, 0, "p", str(pile),\
                self.ylimit, ax)

            figure.text(0.5, 0.04, "base", ha="center")
            figure.text(0.04, 0.5, "coverage", va="center", rotation="vertical")
            matplotlib.pyplot.legend(loc="best")
            matplotlib.pyplot.savefig(self.out_path + str(pile) + ".svg", format='svg', dpi=1200)
            matplotlib.pyplot.close(figure)

#*******************************************************************************

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Plot is a tool for drawing
        different stages of the raven assembler""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("path", help="""Input file in JSON format containing""")
    parser.add_argument("-o", "--out", default=os.getcwd(),
        help="""path in which plotted images will be saved""")

    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument("-y", "--ylimit", help="y axis limit",
        required=True)

    args = parser.parse_args()

    plotter = Plotter(args.path, args.ylimit, args.out)

    with plotter:
        plotter.run()
