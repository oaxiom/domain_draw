domain_draw README

Install:
--------

You will need Python >=3.5.
You will also need matplotlib 

Usage:
------

domain_draw is a tool to draw schematic representations of (typically) protein sequences.

Documentation is a little sparse at the moment.

$ python draw_domains.py --help
Usage: draw_domains.py [options]

Options:
  -h, --help            show this help message and exit
  -i FILENAME           The input filename
  -o OUTPUT_PATH        Path to output the thumbs and full images. Will create
                        own 'thumb' and 'full' directories in whatever you set
                        this to
  -s STYLE, --style=STYLE
                        undocumented style modifier
  -f, --fixed           draw the protein 0 -- 100% (True) or (False) all
                        scaled so that the proteins can be compared in size
  -v, --svg             output 'full' figures as svg files

Input files are expected to be in this format:

>AGAP004733-PA  BTB
    AGAP004733-PA   648     SM00355 481     503     Znf_C2H2        ACCESSORY
    AGAP004733-PA   648     SM00355 564     587     Znf_C2H2        ACCESSORY
    
With a single FASTA entry per protein.

License
-------

Copyright (C) 2011-2020 Andrew Hutchins

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

Except as contained in this notice, the name(s) of the above copyright holders shall not 
be used in advertising or otherwise to promote the sale, use or other dealings in this 
Software without prior written authorization.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
