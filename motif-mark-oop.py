#!/usr/bin/env python
# To be run in motifmark conda environment with only cairo and python installed
import cairo
import argparse
import re

### ARGPARSE SECTION ###
def get_args():
    parser = argparse.ArgumentParser(description="Takes in a fasta file and a text file with one motif per line. Returns an image of the locations of motifs on given sequences.")
    parser.add_argument("-f", "--fasta_file", help="Path to sorted fasta file", type=str, required=True)
    parser.add_argument("-m", "--motifs_file", help="Path to motif file", type=str, required=True)
    parser.add_argument("-t", "--theme", default="vintage", help="Color theme for image. Leave blank for default (vintage). Options are classic, light, dark, mono, contrast, viridi(color blind friendly), and vintage.", type=str, required=False)
    parser.add_argument("-s", "--stagger", default="False", help="Changes formatting of motifs to be below sequence and staggered if True. Default is False", type=str, required=False)
    parser.add_argument("-d", "--dense", default="False", help="Changes formatting of motifs to be solid if True. Default is False. Motifs are more transparent to show overlaps.", type=str, required=False)
    return parser.parse_args()
args = get_args()

### GLOBAL VARIABLES ###
fasta = args.fasta_file
motifsfile = args.motifs_file
DNA = {"A", "T", "C", "G"}
# sets of valid matches to ATCG respectivly
Tmatch = {"U", "W", "K", "Y", "B", "D", "H"}
Amatch = {"W", "M", "R", "D", "H", "V"}
Cmatch = {"S", "M", "Y", "B", "H", "V"}
Gmatch = {"S", "K", "R", "B", "D", "V"}
outputfile = re.findall("([^.]+)", fasta)[0] #matches output file name to input fasta name
outputfile = outputfile + ".png"
### OBJECTS ###
class Sequence:
    '''A genetic sequence'''
    def __init__(self, sequence:str, header: str, length: int):
        ## Data ##
        self.sequence = sequence
        self.header = header
        self.length = length
        self.origin = 0
        self.motiflist = []
        self.motifs = []
        self.exons = []
        self.introns = []
    def __repr__(self):
        rep = 'Sequence for ' + self.header + ''
        return rep
    ## Methods ##
    def place_motifs(self, motifs):
        '''Places motifs from file on sequence, returns as list of tuples (motifname, [location list]).'''
        self.motiflist = motifs
        motiflist = []
        for motif in motifs:
            kmer = len(motif)
            locations = []
            sequence = self.sequence
            for location in range(len(sequence)):
                end = kmer+location
                if end == (len(sequence)+1):
                    break
                curr = sequence[location:end]
                include = translate_motif(motif,curr)
                if include:
                    locations.append(location)
            motiflist.append((motif, locations))
        self.motifs = motiflist

    def place_codingblocks(self):
        '''Assigns list of tuples to exons and introns based on capitalization in sequence. Tuples hold start location and length.'''
        Ilocations = [] #currently introns are not being used, possible future functionality
        Elocations = []
        sequence = self.sequence
        start = 0
        i = 0
        while i < len(sequence):
            if sequence[i].isupper():
                start = i
                while sequence[i].isupper():
                    if i == (len(sequence)-1):
                        break
                    i+=1
                Elocations.append((start, (i-start+1)))
            if sequence[i].islower():
                start = i
                while sequence[i].islower():
                    if i == (len(sequence)-1):
                        break
                    i+=1
                Ilocations.append((start, (i-start+1)))
            i+=1
        self.exons = Elocations
        self.introns = Ilocations

class Theme:
    '''The theme for figure. Includes color scheme, spacing sizes, and transparensy values.'''
    def __init__(self, figsize=int, buffer=int, color=str):
        ## Data ##
        self.figsize = figsize #height of figure
        self.buffer = buffer
        self.colorscheme = color
        self.alpha = 0.8
        self.motifcolors = [self.translate_color([168, 83, 181]), self.translate_color([242, 19, 83]), self.translate_color([103, 168, 45]), self.translate_color([8, 108, 158]), self.translate_color([242, 139, 5])]
        self.backgroundcolor = [1,1,1,0]
        self.strokecolor = [.1, .1, .1,self.alpha]
        self.assign_colors
        self.panelbuffer = 5
        self.stagger = args.stagger
        self.legendsize = 100
        self.itemheight = 25 #height of lables
        self.textbuffer = 7
        self.dense = args.dense

        #self.theme = theme
    def __repr__(self):
        rep = 'Theme:'#include specs here + self.size + 'sequences'
        return rep
    def translate_color(self, color: list): #pycairo takes colors as  0-1 % of rgb values. This translates normal rgb numbers into the correct format.
        newcolor = [color[0]/250, color[1]/250, color[2]/250, self.alpha]
        return newcolor
    def assign_colors(self):
        if self.dense == "True": #makes motifs solid and transparent instead of hollow and more solid.
            self.alpha = 0.6
        if self.colorscheme == "light":
            self.strokecolor = self.translate_color([98, 113, 138])
            self.motifcolors = [self.translate_color([168, 83, 181]), self.translate_color([242, 19, 83]), self.translate_color([103, 168, 45]), self.translate_color([10, 141, 207]), self.translate_color([34, 189, 189])]
        if self.colorscheme == "dark":
            self.backgroundcolor = self.translate_color([14, 23, 36])
            self.strokecolor = self.translate_color([225, 236, 252])
            self.motifcolors = [self.translate_color([168, 83, 181]), self.translate_color([242, 19, 83]), self.translate_color([103, 168, 45]), self.translate_color([10, 141, 207]), self.translate_color([34, 189, 189])]
        if self.colorscheme == "vintage": #modeled after vintage apple colors.
            self.motifcolors = [self.translate_color([97,187, 70]), self.translate_color([253,184,39]), self.translate_color([245,130,31]), self.translate_color([224,58,62]), self.translate_color([150,61,151])]
            self.backgroundcolor = self.translate_color([245,231,206])
            self.strokecolor = self.translate_color([94,41,47])
        if self.colorscheme == "mono": #intended for black and white printing
            self.motifcolors = [self.translate_color([34, 5, 45]), self.translate_color([86, 54, 97]), self.translate_color([119, 92, 127]), self.translate_color([160, 133, 166]), self.translate_color([204, 179, 209])]
        if self.colorscheme == "classic":
            self.motifcolors = [self.translate_color([168, 83, 181]), self.translate_color([242, 19, 83]), self.translate_color([103, 168, 45]), self.translate_color([8, 108, 158]), self.translate_color([242, 139, 5])]
            self.backgroundcolor = [1,1,1,0]
            self.strokecolor = [.1, .1, .1,self.alpha]
        if self.colorscheme == "viridi": #intended to be color blind friendly.
            self.motifcolors = [self.translate_color([162, 12, 193]), self.translate_color([64, 14, 160]), self.translate_color([21, 182, 191]), self.translate_color([60, 226, 10]), self.translate_color([28, 158, 67])]
            self.backgroundcolor = self.translate_color([255, 250, 227])
            self.strokecolor = self.translate_color([105, 81, 56,])
        if self.colorscheme == "contrast": #darker motifs for higher contrast with background and less eye strain.
            self.motifcolors = [self.translate_color([87, 15, 112]), self.translate_color([6, 148, 53]), self.translate_color([212, 135, 11]), self.translate_color([17, 79, 171]), self.translate_color([247, 12, 55])]
            self.backgroundcolor = [1,1,1,0]
            self.strokecolor = [.1, .1, .1,self.alpha]
        self.strokecolor[3] = 0.8

class Figure:
    '''The figure to be drawn.'''
    def __init__(self, sequences:list, theme: Theme):
        ## Data ##
        self.sequences = sequences
        self.theme = theme
        self.longestmotif = 0
        self.legend = (self.longestmotif*self.theme.textbuffer)+50
    def __repr__(self):
        rep = 'Figure for ' + self.size + 'sequences'
        return rep
    ## Methods ##
    def dimensions(self):
        '''Determines sizes for figure based on sizes of sequences. Returns tuple with height of one panel and longest sequence length.'''
        sequences = self.sequences
        buffer = theme.buffer
        labelheight = theme.itemheight+theme.panelbuffer
        longest=0
        # finds longest sequence to scale all sequence graphs
        for item in sequences:
            if longest < len(item.sequence):
                longest = len(item.sequence)
        # There is a lot of math here, it should auto-fit everything but may wig out if buffers or fig size is edited too much.
        # height is the size of the entire sequence figure including label and buffers while panelsize is the section dedicated to the sequence itself.
        panelsize = ((theme.figsize - 2*buffer)/len(sequences)) - labelheight
        height = 2*buffer + (panelsize+labelheight)
        return (height, panelsize, longest)
    def draw_legend(self,length,context):
        '''Creates the color key on the top right of the figure.'''
        num=1
        # Sizes legend
        longestmotif = self.longestmotif
        legend = self.legend
        textbuffer = self.theme.textbuffer
        itemheight = self.theme.itemheight
        motifs = self.sequences[0].motifs
        color = self.theme.strokecolor
        context.set_source_rgba(color[0], color[1], color[2])
        buffer = self.theme.buffer/2
        context.rectangle(length-buffer, buffer, -50-(longestmotif*textbuffer), itemheight*(len(motifs)+1))
        context.stroke()
        # Title
        context.move_to(length-legend/2-buffer-20, buffer+15)
        context.show_text("LEGEND")
        context.stroke()
        # Adds motif names and colors
        for motif in motifs:
            color = self.theme.strokecolor
            context.set_source_rgb(color[0], color[1], color[2])
            context.move_to(length-legend-buffer+5, buffer+itemheight*num+10)
            context.show_text((motif[0]).lower())
            context.stroke()
            color = self.theme.motifcolors[num-1]
            context.set_source_rgb(color[0], color[1], color[2])
            context.rectangle(length-legend-buffer+5+(longestmotif*textbuffer),buffer+itemheight*num+2, legend-10-(longestmotif*textbuffer), itemheight-9)   #legend-(longestmotif*textbuffer)-20,itemheight-9)
            context.fill()
            num+=1
    def draw_figure(self):
        '''Creates final image.'''
        dim = self.dimensions()
        self.longestmotif = len(max(self.sequences[0].motiflist, key = len))
        legend = 50+(self.longestmotif*self.theme.textbuffer)
        self.legend = legend
        itemheight = self.theme.itemheight
        length=(dim[2]+(self.theme.buffer*2))+legend
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, length, theme.figsize)
        context = cairo.Context(surface)
        item = 0 # curr subfigure
        # create top buffer
        color = self.theme.backgroundcolor
        context.set_source_rgb(color[0], color[1], color[2])
        context.rectangle(0, 0, length, theme.buffer)
        context.fill()
        # Fill in sequence section
        for sequence in self.sequences:
            if theme.stagger == "True": #for stagger placement theme
                buffer = theme.buffer
                panel = dim[1]
                sectionbuffer = theme.itemheight+theme.panelbuffer
                section = panel+sectionbuffer
                linelength = len(sequence.sequence)
                linebuffer = theme.buffer + (.5*(dim[2]-linelength)) #used to center lines in page
                center= buffer + item*section + .5*panel-.25*section #denotes veritcal center of curr subfig
                ## create background ##
                color = self.theme.backgroundcolor
                context.set_source_rgb(color[0], color[1], color[2])
                context.rectangle(0, section*item+buffer, length, section)
                context.fill()
                ## create sequence line ##
                context.set_line_width(1)
                color = self.theme.strokecolor
                context.set_source_rgb(color[0], color[1], color[2])
                context.move_to(linebuffer,center)
                context.line_to(linelength+linebuffer,center)
                context.stroke()
                ## Create exons ##
                for exon in sequence.exons:
                    context.set_source_rgba(color[0], color[1], color[2], color[3])
                    context.rectangle(exon[0]+linebuffer,center-10,exon[1],20)
                    context.fill()
                colornum = 0
                ## Create motifs ##
                for motif in sequence.motifs:
                    motif_length = len(motif[0])
                    color = self.theme.motifcolors[colornum]
                    context.set_source_rgba(color[0], color[1], color[2], color[3])
                    for loc in motif[1]:
                        context.rectangle(loc+linebuffer,center+colornum*10+15,motif_length,10)
                        if(theme.dense == "True"):
                            context.fill()
                        else:
                            context.stroke()
                    colornum +=1
            else: #normal placement theme
                buffer = theme.buffer
                panel = dim[1]
                sectionbuffer = theme.itemheight+theme.panelbuffer
                section = panel+sectionbuffer
                linelength = len(sequence.sequence)
                linebuffer = theme.buffer + (.5*(dim[2]-linelength)) #used to center lines in page
                center= buffer + item*section + .5*panel #denotes veritcal center of curr subfig
                ## create background ##
                color = self.theme.backgroundcolor
                context.set_source_rgb(color[0], color[1], color[2])
                context.rectangle(0, section*item+buffer, length, section)
                context.fill()
                ## create sequence line ##
                context.set_line_width(1)
                color = self.theme.strokecolor
                context.set_source_rgb(color[0], color[1], color[2])
                context.move_to(linebuffer,center)
                context.line_to(linelength+linebuffer,center)
                context.stroke()
                ## Create exons ##
                for exon in sequence.exons:
                    context.set_source_rgba(color[0], color[1], color[2], color[3])
                    context.rectangle(exon[0]+linebuffer,center-20,exon[1],40)
                    context.fill()
                colornum = 0
                ## Create motifs ##
                for motif in sequence.motifs:
                    motif_length = len(motif[0])
                    color = self.theme.motifcolors[colornum]
                    context.set_source_rgba(color[0], color[1], color[2], color[3])
                    for loc in motif[1]:
                        context.rectangle(loc+linebuffer,center-30+colornum,motif_length,60-colornum*2)
                        if(theme.dense == "True"):
                            context.fill()
                        else:
                            context.stroke()
                    colornum +=1
            ## Create labels ##
            color = self.theme.strokecolor
            context.set_source_rgba(color[0], color[1], color[2])
            context.move_to(buffer, buffer + panel*(item+1) + sectionbuffer*(item+0.5))
            context.show_text(sequence.header)
            context.stroke()
            item +=1
        self.draw_legend(length, context)
        ## create bottom buffer ##
        color = self.theme.backgroundcolor
        context.set_source_rgb(color[0], color[1], color[2])
        context.rectangle(0, theme.buffer+(dim[1]+itemheight+theme.panelbuffer)*len(self.sequences), length, theme.buffer)
        context.fill()
        surface.write_to_png(outputfile)
        surface.finish()

### GLOBAL FUNCTIONS ###
def import_fasta_file(filepath: str):
    '''opens file from string path, returns list of strings of sequences contained. First item in list is the number of sequences contained.'''
    with open(filepath, "r") as file:
        output=[]
        curr = 0
        for line in file:
            line = line.strip()
            if line.startswith(">"): #if a header
                new = [line, ""]
                output.append(new)
                curr+=1
            else:
                output[curr-1][1] += line
        return(output)

def import_motif_file(filepath: str):
    '''opens file from string path, returns list of strings of motifs contained.'''
    with open(filepath, "r") as file:
        output=[]
        for line in file:
            line = line.strip()
            output.append(line)
        return(output)

def create_sequences(filelist: list):
    '''takes in output from import_fafsa_file and returns same list with sequence strings replaced with sequence objects'''
    num = len(filelist)
    for i in range(0, num):
        newobject = (Sequence(filelist[i][1], filelist[i][0], len(filelist[i][1])))
        filelist[i] = newobject
    return(filelist)

def translate_motif(motif, section):
    '''Check if each letter of give kmer is acceptable for the motif. Section and motif must be same length. Return true or false'''
    # Capitalize everyhting
    motif = motif.upper()
    section = section.upper()
    motifs = []
    # See if motif letter is one of the matches for given letter, based on the sets in the global variables section.
    for i in range(len(section)):
        letter = section[i]
        curr = motif[i]
        if not letter in DNA:
            return(False)
        if(letter == curr or curr == "N"):
            continue
        else:
            if letter == "T":
                if not curr in Tmatch:
                    return(False)
                else:
                    continue
            if letter == "A":
                if not curr in Amatch:
                    return(False)
                else:
                    continue
            if letter == "C":
                if not curr in Cmatch:
                    return(False)
                else:
                    continue            
            if letter == "G":
                if not curr in Gmatch:
                    return(False)
                else:
                    continue
    return(True)


## Run the programs ##
motifs = import_motif_file(motifsfile)
output = create_sequences(import_fasta_file(fasta))
for item in output:
    item.place_motifs(motifs)
    item.place_codingblocks()
theme = Theme((len(output)*100+100), 25, args.theme) #figure size can be changed but this autofits nicely.
theme.assign_colors()
fig = Figure(output, theme)
fig.draw_figure()