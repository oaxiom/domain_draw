"""

draw domain alignments

"""

import sys, os, gc, random
from optparse import OptionParser
import gzip as gzipfile

from .data import *
from .tools import *
from .adjustText import adjust_text

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

class schematic:
    def __init__(self, pdf=True, fixed=True):
        """
        **Purpose**
            initiator

        """
        self.pdf = pdf
        self.fixed = fixed
        self.col_map = {}

    def parse_file(self, domain_filename,
        gzip=False,
        fixed=True, pdf=False, output_path=None, **kargs):
        """
        **Purpose**
            Entry point for file parsing

        **Arguments**
            domain_filename
                filename to parse.


            format (Required),
                Format: 'dfam'



                Format: 'fasta_style'
                Expected format:
                >AGAP004733-PA  BTB
                    AGAP004733-PA   648     SSF54695        1       144     SSF54695        FAMILY-DEFINING
                    AGAP004733-PA   648     SM00355 481     503     Znf_C2H2        ACCESSORY
                    AGAP004733-PA   648     SM00355 564     587     Znf_C2H2        ACCESSORY

                ><name> <family_type>
                    [List of domains]
                    <id>    <length of protein> <HMM ID>    <left>  <right> <name>  <type>

            gzip (Optional, default=False)
                Is the domain_filename gzipped?

            fixed (Optional, default=True)
                draw all schemas with a fixed length.
                IF set to False, the longest peptide will be scaled to the full widht
                and everything else drawn in relation to that.
                This is not guaranteed to give sensible results.

            pdf (Optional, default=False)
                Save image as a pdf

            output_path (Optional, default=None)
                path to export the files to.
                Will be created if not present (Not guaranteed to work).
        """
        valid_file_formats = set(['fasta_style', 'dfam'])
        assert format in valid_file_formats, '{0} not in {1} valid_file_formats'.format(format, valid_file_formats)

        if output_path:
            if not os.path.exists(output_path):
                os.mkdir(output_path)
            self.output_path = output_path

        self.fixed = fixed
        self.pdf = pdf

        self.data = []

        if format == 'fasta_style':
            self.data = self.__load_fasta_style(filename, gzip)
        elif format == 'dfam':
            self.data = self.__load_dfam_style(filename, gzip)

    def __load_fasta_style(self, filename, gzip=False):
        """


        """
        pass

    def __load_fasta_style(self, filename, gzip=False):
        '''
                        Expected format:
                >AGAP004733-PA  BTB
                    AGAP004733-PA   648     SSF54695        1       144     SSF54695        FAMILY-DEFINING
                    AGAP004733-PA   648     SM00355 481     503     Znf_C2H2        ACCESSORY
                    AGAP004733-PA   648     SM00355 564     587     Znf_C2H2        ACCESSORY

                ><name> <family_type>
                    [List of domains]
                    <id>    <length of protein> <HMM ID>    <left>  <right> <name>  <type>
        '''
        oh = open(filename, "rt")
        # files apear to be a bit like FASTA files, then tab separated elements

        item = None
        types = []

        self.max_len = 0 # get the maximum length of the peptides
        for line in oh:
            if line:
                if ">" in line:
                    # append the last item to the list
                    if item:
                        self.data.append(item)
                    # Get a new item
                    try:
                        type = " ".join(line.split()[1:])
                        if "/" in type:
                            type = "Mixed:%s" % type
                    except IndexError: # There is no type
                        type = None

                    item = {"name": line.strip().split()[0].replace(">", ""),
                        "type": type, "domains": [], "len": 1}
                else:
                    l = line.split()
                    if len(l) == 7:
                        item["domains"].append({"name": l[5], "pos": (int(l[3])-1, int(l[4])-1), "db": l[2], "fam": l[6]}) # The display is zero-ordered. Correct positions.
                    elif len(l) == 6: # the family/specific column is missing. Seen in the wild occasionally
                        item["domains"].append({"name": l[5], "pos": (int(l[3])-1, int(l[4])-1), "db": l[2]}) # The display is zero-ordered. Correct positions.

                    item["len"] = int(l[1])
                    if self.max_len < item["len"]:
                        self.max_len = item["len"]

        if item: # Make sure the last item gets added
            self.data.append(item)

        oh.close()

    def draw_all(self, style="ubl", thumbs=False):
        """
        **Purpose**
            Just a simple helper function when you want to draw all motifs from the file

        **Arguments**
            style
                "ubl" - the UBL domain E2,E3 style
                "pfsmsff" - PFAM, SFF, SMART colours
        """
        if style == "ubl":
            self.col_map = ubl_col_map
        elif style == "pfsmsff":
            self.col_map = acast_col_map
        elif style == "unk_domains":
            self.col_map = unk_domains
        elif style == "ptp":
            self.col_map = ptp_map

        for n, item in enumerate(self.data):
            if style == "ubl":
                self.draw(item, self.__draw_ubl_style, thumbs)
            elif style == "pfsmsff":
                self.draw(item, self.__draw_db_style, thumbs)
            elif style == "gen":
                self.draw(item, self.__draw_gen_style, thumbs)
            elif style == "unk_domains":
                self.draw(item, self.__draw_unk_style, thumbs)
            elif style == "ptp":
                self.draw(item, self.__draw_ubl_style, thumbs)

        return(None)

    def draw(self, item, filename, style, thumb=False):
        """
        **Purpose**
            draw the item.

        **Arguments**
            item
                A dictionary in the form:
                {"name": <name>,
                "type": <type>,
                "domains: [<list of domains> {"name": <name>, "pos": (l, r), "fam": <family>},
                "len": <length of protein>
                }

        **Returns**
            None and a file in <item>.png
        """
        valid_styles = {
            #'dudedb': self.__draw_ubl_style, # need revision
            #'ptp': self.__draw_ubl_style, # need revision
            'episcan': self.__draw_gen_style,
            }

        assert style in valid_styles, '{0} style not in valid_styles ({1})'.format(style, valid_styles)

        if thumb:
            p = {"pad1": item["len"] * 0.02,
                "pad2": item["len"] * 0.003,
                "figsize": (1.5, 0.3),
                "lpad": 0.16,
                "evpad": 0.5,
                "nvpad": 0.4}
        else:
            if not self.fixed:
                p = {"pad1": self.max_len * 0.035, # pad out 1% of the figure left and right
                    "pad2": self.max_len * 0.003,
                    "figsize": (25, 1),
                    "lpad": 0.02, # vertical size of the central line
                    "evpad": 0.25, # vertical size of the 'normal' rectangle domain boxes.
                    "nvpad": 0.2,
                    "titlepos": "left",
                    "titlesize": 13} # vertical size of the 'enhanced' rectangle domain boxes.
            else: # ie fixed
                p = {"pad1":  item["len"] * 0.03, # pad out 1% of the figure left and right
                    "pad2": item["len"] * 0.003,
                    "figsize": (7, 1),
                    "lpad": 0.02, # vertical size of the central line
                    "evpad": 0.20, # vertical size of the 'normal' rectangle domain boxes.
                    "nvpad": 0.2,
                    "titlesize": 13
                    } # vertical size of the 'enhanced' rectangle domain boxes.


        fig = plt.figure(figsize=p["figsize"])

        ax = fig.add_subplot(211)
        axlab = fig.add_subplot(212)

        ax.set_position([0, 0.5, 1.0, 0.5])
        axlab.set_position([0, 0.0, 1.0, 0.60])

        ax.add_patch(Rectangle((0,-p["lpad"]), item["len"]-1, p["lpad"]*2, ec="none", fc="black", color="black")) # The line for the protein

        if self.fixed:
            ax.set_xlim([-p["pad1"], item["len"]+p["pad1"]])
            axlab.set_xlim([-p["pad1"], item["len"]+p["pad1"]])
        else:
            ax.set_xlim([-p["pad1"], self.max_len+p["pad1"]])
            axlab.set_xlim([-p["pad1"], item["len"]+p["pad1"]])

        ax.set_ylim([-1, 1])

        axlab.set_ylim([0, 10])
        axlab.set_facecolor('none')

        if len(item["domains"]) > 0:
            valid_styles[style](ax, axlab, item, p, thumb)

        if not thumb: # numbers showing the aa position of peptide
            tpad = p['pad1']/10
            ax.text(-tpad, 0, str(0), ha="right", va="center", fontsize=5, color="black", zorder=100001)
            ax.text(item["len"]+tpad, 0, item['len'], ha="left", va="center", fontsize=5, color="black", zorder=100001)

        ax.set_xticklabels("")
        ax.set_yticklabels("")
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.tick_params(left=None, top=None, bottom=None, right=None)
        ax.set_frame_on(False)

        axlab.set_xticklabels("")
        axlab.set_yticklabels("")
        axlab.set_ylabel("")
        axlab.set_xlabel("")
        axlab.tick_params(left=None, top=None, bottom=None, right=None)
        axlab.set_frame_on(False)

        fig.savefig(filename)
        plt.close(fig) # Free up the memory

    def __draw_ubl_style(self, ax, item, p, thumb=False):
        """
        drawing style for the ubiquitin ligase database
        """
        low_labs = []
        for d in item["domains"]:
            if d["fam"] == "FAMILY-DEFINING":
                # special code to convert RING finger -> RNF finger
                item["type"] = item["type"].replace("RING finger", "RNF finger")

                if item["type"] in self.col_map:
                    col = self.col_map[item["type"]]
                    label = item["type"]
                elif d["db"] in domain_label_lookup and domain_label_lookup[d["db"]] in self.col_map: # This is where mixed domains will end up.
                    col = self.col_map[domain_label_lookup[d["db"]]]
                    label = domain_label_lookup[d["db"]]
                elif d["name"] in self.col_map:
                    col = self.col_map[d["name"]]
                    label = d["name"]
                else:
                    print("Warning: '%s' not found in colour map" % d["db"])
                    col = "pink"
                    label = d["db"]

                ax.add_patch(Rectangle((d["pos"][0], -p["evpad"]), d["pos"][1] - d["pos"][0], p["evpad"]*2,
                    ec="none", fc=col, lw=0.5, zorder=100000))

                if not thumb: # numbers showing the aa position of the domain
                    ax.text(d["pos"][0]+p["pad2"], 0, str(d["pos"][0]+1), ha="left", va="center", fontsize=5, color="black", zorder=100001)
                    ax.text(d["pos"][1]-p["pad2"], 0, str(d["pos"][1]+1), ha="right", va="center", fontsize=5, color="black", zorder=100001)
            else:
                ax.add_patch(Rectangle((d["pos"][0], -p["nvpad"]), d["pos"][1] - d["pos"][0], p["nvpad"]*2,
                    ec="none", fc="grey", lw=0.5, zorder=2))

            if not thumb:
                if d["fam"] == "FAMILY-DEFINING":
                    l = {"p": (d["pos"][0] + d["pos"][1])/2, "lab": str(label), "fs": 9, "col": "black"}
                else:
                    l = {"p": (d["pos"][0] + d["pos"][1])/2, "lab": str(d["name"]), "fs": 6, "col": "grey"}
                low_labs.append(l)

        # sort out colliding labels:
        for l in low_labs:
            ax.text(l["p"], -0.5, l["lab"], ha="center", va="center", fontsize=l["fs"], color=l["col"])

        if not thumb: # the title
            ax.text(item["len"]/2, 0.7, "%s (%s)" % (item["name"], item["type"]), color="black", fontsize=p["titlesize"], ha="center")

    def __draw_db_style(self, ax, item, p, thumb=False):
        """
        a more generic style of drawing. This one colours the domain depending upon which
        db it came from
        """
        for d in item["domains"]:
            if "SSF" in d["db"]:
                type = "SUPERFAMILY"
            elif "PF" in d["db"]:
                type = "Pfam-A"
            elif "SM" in d["db"]:
                type = "SMART"
            else:
                print("Warning: '%s' not found in colour map" % d["db"])

            ax.add_patch(Rectangle((d["pos"][0], -0.25), d["pos"][1] - d["pos"][0], 0.5,
                ec="black", fc=self.col_map[type], lw=0.5))
            t = ax.text((d["pos"][0] + d["pos"][1])/2, -0.5, str(d["name"]), ha="center", va="center", fontsize=6, color="black")
            texts.append(t)

        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

        ax.text(item["len"]/2, 0.7, "ID: {0}".format(item["name"]), color="black", fontsize=p["titlesize"], ha="center")

    def __draw_gen_style(self, ax, axlab, item, p, thumb=False):
        """
        The most generic drawing style.
        This one colours the domain depending upon which
        db it came from
        """
        texts = []

        for d in item["domains"]:
            ax.add_patch(Rectangle((
                d["pos"][0], -0.25), d["pos"][1] - d["pos"][0], 0.5,
                ec="lightgrey",
                fc='grey', # self.col_map[d["name"]],
                lw=0.5))

            if not thumb:
                r = random.randint(0, 10) / 10 # Give adjust_text something to work with
                t = axlab.text((d["pos"][0] + d["pos"][1])/2, 9-r, str(d["name"]), ha="center", va="center", fontsize=6, color="black")
                texts.append(t)

            if not thumb: # numbers showing the aa position of the domain
                ax.text(d["pos"][0]+p["pad2"], 0.4, str(d["pos"][0]+1), ha="center", va="center", fontsize=5, color="black", zorder=100001)
                ax.text(d["pos"][1]-p["pad2"], 0.4, str(d["pos"][1]+1), ha="center", va="center", fontsize=5, color="black", zorder=100001)

        adjust_text(texts,
            ax=axlab,
            #arrowprops=dict(arrowstyle="-", color='k', lw=0.5),
            only_move={'text': 'y'},
            autoalign='y',
            text_from_points=False,
            force_text=0.2,
            ha="center",
            va="center"
            )

    def __draw_unk_style(self, ax, item, p, thumb=False):
        """
        This is the drawing style for the domains in the CD4+
        T cell paper figure.
        """
        a = True
        for i, d in enumerate(item["domains"]):
            if d["name"] not in self.col_map:
                print("Warning: '%s' not found in colour map" % d["name"])
            else:
                ax.add_patch(Rectangle((d["pos"][0], -0.25), d["pos"][1] - d["pos"][0], 0.5,
                    ec="none", fc=self.col_map[d["name"]], lw=0.5))

            if a:
                ax.text((d["pos"][0] + d["pos"][1])/2, -0.5, str(d["name"]), ha="center", va="center", fontsize=6, color="black")
            else:
                ax.text((d["pos"][0] + d["pos"][1])/2, -0.8, str(d["name"]), ha="center", va="center", fontsize=6, color="black")
            a = not a

            #ax.text((d["pos"][0] + d["pos"][1])/2, -0.7, str(d["db"]), ha="center", va="center", fontsize=6, color="black")

        ax.text(0, 0.5, "%s (%s amino acids)" % (item["name"],  item["len"]), color="black", fontsize=20, ha="left")

if __name__ == "__main__":
    parser = OptionParser() # Should be replaced with argparse
    parser.add_option("-i",
        dest="filename", action="store",
        help="The input filename")
    parser.add_option("-o",
        dest="output_path", default=".",
        help="Path to output the thumbs and full images. Will create own 'thumb' and 'full' directories in whatever you set this to")
    parser.add_option("-s", "--style",
        dest="style", default="ubl",
        help="undocumented style modifier")
    parser.add_option("-f", "--fixed",
        dest="fixed", action="store_true", default=False,
        help="draw the protein 0 -- 100% (True) or (False) all scaled so that the proteins can be compared in size")
    parser.add_option("-v", "--svg",
        dest="svg", action="store_true", default=False,
        help="output 'full' figures as svg files")
    parser.add_option("-c", "--collate", default=False,
        dest="collate", action="store_true",
        help="Scan the domain data and list the FAMILY-DEFINING domains")

    (options, args) = parser.parse_args()

    if options.collate:
        # scan the file and collect all of the FAMILY-DEFINING categories.
        collate_family_defining(options.filename)
    else:
        t = schematic(options.filename, output_path=os.path.join(options.output_path, "full"), fixed=options.fixed, svg=options.svg)
        t.draw_all(style=options.style, thumbs=False)

        t = schematic(options.filename, output_path=os.path.join(options.output_path, "thumbs"))
        t.draw_all(style=options.style, thumbs=True)
