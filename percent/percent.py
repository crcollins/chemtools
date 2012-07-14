import argparse
import sys

parser = argparse.ArgumentParser(description="This program generates homo/lumo fits and graphs.")
parser.add_argument('files', metavar='file', type=str, nargs='*', help='The name of single file.')
parser.add_argument('-o', metavar='output', action="store", dest="output", type=str, help='The name of the file to ouput to.')
parser.add_argument('-d', action="store", dest="draw", default='', help='Toggles drawing out the selected points.')

parser.add_argument('-C', action="store_true", dest="center", default=False, help='Toggles adding the center to total.')
parser.add_argument('-V', action="store_true", dest="vertical", default=False, help='Toggles including the center in the vertical.')
parser.add_argument('-H', action="store_true", dest="horizontal", default=False, help='Toggles including the center in the horizontal.')
parser.add_argument('-Z', action="store_true", dest="zero", default=False, help='Toggles excluding the zero position value.')
parser.add_argument('-E', action="store_true", dest="error", default=False, help='Toggles showing error messages.')

class Output(object):
    def __init__(self, args):
        self.args = args
        self.errors = []

    def write(self, line):
        if self.outputfile:
            self.outputfile.write(line)
        else:
            print line.rstrip('\n')

    def write_file(self):
        if self.args.output:
            self.outputfile = open(self.args.output, 'w')
        else:
            self.outputfile = None

        for filename in self.args.files:
            try:
                self.write(filename + "\n")
                self.run(filename)
                self.write('\n')
            except Exception as ex:
                self.errors.append(ex)

        if self.args.error:
            self.write("\n---- Errors (%i) ----" % len(self.errors))
            for error in self.errors:
                self.write(str(error) + "\n")
        if self.outputfile:
            self.outputfile.close()
        else:
            raw_input("<Press Enter>")

    def draw(self, xvals, yvals, vert, horz, center):
        try:
            import Image
            import ImageDraw
        except ImportError:
            return
        xmin, xmax = min(xvals), max(xvals)-min(xvals)
        ymin, ymax = min(yvals), max(yvals)-min(yvals)
        xvals = [(x-xmin)/xmax for x in xvals]
        yvals = [(y-ymin)/ymax for y in yvals]

        xres = 400
        yres = int(xres * xmax / ymax)
        img = Image.new("RGB", (xres+1, yres+1))
        draw = ImageDraw.Draw(img)

        for i, (x, y) in enumerate(zip(xvals, yvals)):
            # center the values around zero
            xnorm = x - .5
            ynorm = y - .5

            if i in center:
                color = (255, 0, 0)
            elif i in vert:
                color = (0, 255, 0)
            elif i in horz:
                color = (150, 170, 255)
            else:
                color = (255, 255, 255)

            draw.ellipse((x * xres - 4, y * yres - 4, x * xres + 4, y * yres + 4), fill=color)
        img.save(self.args.draw)
        img.show()

    def split_groups(self, xvals, yvals):
        '''Returns the index values for each of the groups '''
        xmin, xmax = min(xvals), max(xvals)-min(xvals)
        ymin, ymax = min(yvals), max(yvals)-min(yvals)
        xvals = [(x-xmin)/xmax for x in xvals]
        yvals = [(y-ymin)/ymax for y in yvals]

        center = []
        vert = []
        horz = []
        radiussq = (1 / 8.0) ** 2
        for i, (x, y) in enumerate(zip(xvals, yvals)):
            xnorm = x - .5
            ynorm = y - .5
            if xnorm ** 2 + ynorm ** 2 < radiussq:
                center.append(i)
            elif ynorm > abs(xnorm) or ynorm < -abs(xnorm):
                horz.append(i)
            elif xnorm > abs(ynorm) or xnorm < -abs(ynorm):
                vert.append(i)
            else:
                raise ValueError
        return vert, horz, center

    def get_coords(self, filename):
        '''Extracts all the coords from the log file. '''
        nums = []
        xvals = []
        yvals = []
        with open(filename, 'r') as f:
            state = 0
            for line in f:
                if state == 0 and "Number     Number       Type             X           Y           Z" in line:
                    state = 1
                elif state == 1 and "-"*60 in line:
                    state = 2
                elif state == 2:
                    if "-"*60 in line:
                        state = 3
                    else:
                        s = [float(x) for x in line.split()[3:5]]
                        xvals.append(s[0])
                        yvals.append(s[1])
        return xvals, yvals

    def get_numbers(self, filename):
        with open(filename, 'r') as f:
            state = 0
            homo = [[]]
            lumo = [[]]
            for line in f:
                split = line.split()
                if state == 0 and "O         O         O         O         O" in line:
                    state = 1
                elif state == 1:
                    if len(split) == 9:
                        homo.append([])
                    if "V         V         V         V         V" in line:
                        state = 2
                    elif len(split) != 5:
                        homo[-1].append(float(split[-1]))
                elif state == 2:
                    if len(split) == 9:
                        lumo.append([])
                    if "Density Matrix:" in line:
                        state = 3
                    else:
                        lumo[-1].append(float(split[-5]))
        return homo, lumo

    def calc_groups(self, value, vert, horz, center):
        # drop the first because wrong?
        valuesq = [sum(y ** 2 for y in x) for x in value[1:]]
        if not self.args.zero:
            valuesq[0] += value[0][0]**2
        if self.args.vertical:
            vert += center
        if self.args.horizontal:
            horz += center

        vertvalue = sum(x for i, x in enumerate(valuesq) if i in vert)
        horzvalue = sum(x for i, x in enumerate(valuesq) if i in horz)
        centervalue = sum(x for i, x in enumerate(valuesq) if i in center)

        totalvalue = sum(valuesq)
        if not self.args.center:
            totalvalue -= centervalue

        return totalvalue, vertvalue, horzvalue, centervalue

    def disp(self, name, total, vert, horz, center):
        self.write('%s              Raw     Percent\n' % name)
        self.write("Vertical:    %f   %f\n" % (vert, (vert / total * 100)))
        self.write("Horizontal:  %f   %f\n" % (horz, (horz / total * 100)))
        if self.args.center:
            self.write("Center:      %f   %f\n" % (center, (center / total * 100)))
        self.write("Total:       %f\n\n" % (total))

    def run(self, filename):
        xvals, yvals = self.get_coords(filename)
        vert, horz, center = self.split_groups(xvals, yvals)
        if self.args.draw:
            self.draw(xvals, yvals, vert, horz, center)
        homo, lumo = self.get_numbers(filename)

        self.disp("HOMO", *self.calc_groups(homo, vert, horz, center))
        self.disp("LUMO", *self.calc_groups(lumo, vert, horz, center))


if __name__ == '__main__':
    if len(sys.argv) > 1:
        a = Output(parser.parse_args(sys.argv[1:]))
        a.write_file()
    else:
        args = raw_input('Arguments: ')
        a = Output(parser.parse_args(args.strip().split()))
        a.write_file()
