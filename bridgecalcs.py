#INPUT SPECS: nested list with each individual list representing one shape
#Each shape should have the following information:
#[b (x direction dimension), h (y direction dimension), x position of bottom left corner, y position of bottom left corner]
#x, y axis goes from 0 - 200 mm

#For horizontal shapes above the centroid, add a list of the buckling cases it faces
#e.g. for _____________________ it would be [2, 1, 1, 2]
#           |      |      |
#for ________________ it would be [2]
#    |              | 

shapes = [
    #bottom flange
    [80, 1.27, 20, 0], 
    #webs
    [1.27, 20+96.19-1.27, 20, 1.27], 
    [1.27, 20+96.19-1.27, 100-1.27, 1.27],

    #glue joints
    [10, 1.27, 20, 20+96.19],
    [10, 1.27, 90, 20+96.19],
    

    #top flange
    [112, 2.54, 4, 20+96.19 + 1./27, [2, 1, 2]]

]

#GLUE INFO: [Height of   glue (same axis as shape), width of glued section]
glue = [
    [96.19 + 1.27 * 2 + 20, 6.27*2],
    [96.19 + 1.27 + 20, 6.27*2],
]

#NUMBER OF DIAPHRAGMS 
diaphragms = 8


#Length of bridge - max possible length -> biggest seperation between diaphragms -> more conservative calc of FOS
L = 1270

#Max moment (N mm):
#M = 23980
M = 78099
#M =  69445.33333333


#Max shear (N)
V = 284.517
#V = 240

mu = 0.2

E = 4000


#==================================#
#DON'T CHANGE ANYTHING BELOW THIS!!!
#==================================#


import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

def draw_cross_section(shapes):
    fig, ax = plt.subplots(figsize=(8, 8))
    
    def save_on_close(event):
        fig.savefig("crosssec.png", format="png")

    fig.canvas.mpl_connect('close_event', save_on_close)

    # Set axis limits
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 200)
    
    # draw centroid
    for shape in shapes:
        width, height, x, y, a = shape

        # plot rectangle
        rect = patches.Rectangle(
            (x, y), width, height, linewidth=1, edgecolor='blue', facecolor='lightblue', alpha=0.7
        )
        ax.add_patch(rect)
        
        # label length of each rectangle at its center
        x_centroid, y_centroid = localcentroid(shape)
        ax.text(x_centroid, y_centroid, f"{max(width, height):.2f}", 
                color="red", fontsize=8, ha="center", va="center")
    
    ax.set_title("Cross Section Visual")
    ax.set_xlabel("X Position (mm)")
    ax.set_ylabel("Y Position (mm)")
    ax.grid(True)
    
    #plot centoridal axis
    ybar = centroid(shapes)
    ax.axhline(y=ybar, color='blue', linestyle='--', linewidth=1, label=f"Centroidal Axis (y ~ {ybar:.3f})")

    for height, width in glue:
        ax.axhline(y=height, color='green', linestyle='--', linewidth=1, label=f"Glued Surface (Width = {width:.2f}, y = {height:.2f})")

    plt.legend(loc = "upper left")

    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()



def convert(shapes):
    for i in range(len(shapes)):
        if len(shapes[i]) != 5:
            shapes[i].append([])
    return shapes


#Converts from orig data structure to 
#[b, h, x, y, [cases]]
#(old), data structure changed
def convertOld(shapes):
    for i in range(len(shapes)):
        if shapes[i][0] == 0:
            shapes[i][0] = shapes[i][1]
            shapes[i][1] = 1.27
        elif shapes[i][0] == 90:
            shapes[i][0] = 1.27
        else:
            raise ValueError("Invalid rotation: must be 0 or 90.")
        if len(shapes[i]) != 5:
            shapes[i].append([])
    return shapes

#helper func finds local cross sec of rectangles
def localcentroid(shape):
    width, height, x, y, a = shape
    return x + width/2, y + height/2

def centroid(shapes):
    ybar = 0
    totalarea = 0
    for shape in shapes:
        x, y = localcentroid(shape)
        area = shape[0] * shape[1]
        ybar += area * y
        totalarea += area
    ybar /= totalarea
    return ybar

def secondmomentarea(shapes, ybar):
    total = 0
    for shape in shapes:
        x, y = localcentroid(shape)
        total += (shape[0] * shape[1] ** 3)/12
        total += (shape[0] * shape[1]) * (y - ybar) ** 2
    return total

#splits cross sec into shapes above a given line
def splitAbove(shapes, line):
    newshapes = []
    for shape in shapes:
        b, h, x, y, a = shape
        if y + h < line:
            continue
        if y < line:
            newshapes.append([b ,(y + h - line), x, line, []])
        else:
            newshapes.append(shape)
    return newshapes

#helper func returns total area of shapes
def area(shapes):
    area = 0
    for shape in shapes:
        area += shape[0] * shape[1]
    return area

#computes flexural comp of cross sec
def flexuralCompression(shapes):
    ybar = centroid(shapes)
    top = 0
    #find max h + y
    for shape in shapes:
        if shape[3] + shape[1] > top:
            top = shape[3] + shape[1]
    return (M * (top - ybar))/secondmomentarea(shapes, ybar)

#computes flexural tens of cross sec
def flexuralTension(shapes):
    ybar = centroid(shapes)
    bottom = 200
    #find max h + y
    for shape in shapes:
        if shape[3] < bottom:
            bottom = shape[3]
    return (M * (ybar - bottom))/secondmomentarea(shapes, ybar)

#computes glue shear
#t = VQ/Ib
def glueShear(shapes, glue):
    ybar = centroid(shapes)
    I = secondmomentarea(shapes, ybar)
    taumax = 0
    for y, b in glue:
        glueshapes = splitAbove(shapes, y)
        Q = area(glueshapes) * (centroid(glueshapes) - ybar)
        if taumax < abs(V*Q/(I*b)):
            taumax = abs(V*Q/(I*b))
    return taumax

#computes shear at centroid
def centroidShear(shapes):
    ybar = centroid(shapes)
    I = secondmomentarea(shapes, ybar)
    areaAboveCentroid = splitAbove(shapes, ybar)
    Q = area(areaAboveCentroid) * (centroid(areaAboveCentroid) - ybar)
    b = 0
    for shape in areaAboveCentroid:
        if shape[3] == ybar:
            b += shape[0]
    return V*Q/(I*b)

#helper func splits horizontal section into multiple smaller ones
def splitVertical(shape, xvals):
    if not xvals:
        return []
    b,h,x,y,a = shape
    if xvals[0] != x:
        xvals.insert(0, x)
    if xvals[-1] != x+b:
        xvals.append(x+b)
    ans = []
    for i in range(1, len(xvals)):
        ans.append([xvals[i] - xvals[i-1], h, y + h])
    return ans

#helper func to split cross sec into lists of case1, case2, and case3 plate buckling
def splitCases(shapes):
    ybar = centroid(shapes)
    shapes = splitAbove(shapes, ybar)

    case1 = []
    case2 = []
    case3 = []

    for shape in shapes:
        #If section is vertical -> case 3
        if shape[0] == 1.27:
            b, h, x, y, a = shape
            case3.append([b, h, y + h])
            continue
        #if list of cases is not empty:
        if shape[4]:
            #find vertical pieces closest to the horizontal one and store them
            dist = 10000
            verticals = []
            for i in shapes:
                if i[0] == 1.27 and shape[3] - (i[3] + i[1]) < dist:
                    dist = shape[3] - (i[3] + i[1])
                    verticals = []
                    verticals.append(i)
                elif i[0] == 1.27 and shape[3] - (i[3] + i[1]) == dist:
                    verticals.append(i)
            #split the shape into new shapes according to its cases, locations of verticals
            xvals = []
            for i in verticals:
                x, y = localcentroid(i)
                xvals.append(x)
            cases = splitVertical(shape, xvals)
            for i in range(len(shape[4])):
                if shape[4][i] == 1:
                    case1.append(cases[i])
                else:
                    case2.append(cases[i])
    return case1, case2, case3

#computes buckling stress of case1 shapes, returns minimum
def case1(shapes, I, ybar):
    stresses = []
    fos = []
    for b, t, y in shapes:
        tau = (4*(math.pi**2)*E) / (12*(1-mu**2)) * (t / b) ** 2
        stress = (M * (y-ybar))/I
        stresses.append(tau)
        fos.append(tau/stress)
    return stresses, fos

#computes buckling stress of case2 shapes, returns minimum
def case2(shapes, I, ybar):
    stresses = []
    fos = []
    for b, t, y in shapes:
        tau = (0.425*(math.pi**2)*E) / (12*(1-mu**2)) * (t / b) ** 2
        stress = (M * (y-ybar))/I
        stresses.append(tau)
        fos.append(tau/stress)
    return stresses, fos


#computes buckling stress of case3 shapes, returns minimum
#****ASSUMES THERE IS ALWAYS A GLUE TAB ON TOP OF VERTICAL WEB
#FOR NO GLUE TAB REMOVE THE + 1.27 IN EQUATION
def case3(shapes, I, ybar):
    stresses = []
    fos = []
    for t, b, y in shapes:
        tau = (6*(math.pi**2)*E) / (12*(1-mu**2)) * (t / (b + 1.27)) ** 2
        stress = (M * (y - ybar + 1.27))/I
        stresses.append(tau)
        fos.append(tau/stress)
    return stresses, fos


#finds max height of the cross section 
#helpfer func for case4 buckling
def height(shapes):
    h_min = 100000
    h_max = 0
    for b, h, x, y, a in shapes:
        if y < h_min:
            h_min = y
        if y+h > h_max:
            h_max = y + h
    return h_max - h_min

#finds shear buckling stress
def case4(shapes):
    t = 1.27
    h = height(shapes)
    a = L / (diaphragms + 1)
    tau = (5*(math.pi**2)*E) / (12*(1-mu**2)) * ((t / h) ** 2 + (t/a) **2)
    return tau

#outputs everything beautifully and nicely tabulated
#adjust spacing var to change space between label and value
#!!Aditionally performs some minor calculations i.e. determining FOS'
#!! takes the minimum between 4 and shear buckling stress calculated in case4 as the allowable shear stress
def bundledOutput(shapes):
    shapes = convert(shapes)
    spacing = 50
    ybar = centroid(shapes)
    I = secondmomentarea(shapes, ybar)

    print("\n" + "-" * 80)
    print(f"{'Centroid:':<{spacing}} {ybar}")
    print(f"{'I:':<{spacing}} {I}")
    print("-" * 80)

    print("\n" + "-" * 80)
    print(f"{'Flexural Tension:':<{spacing}} {flexuralTension(shapes)}")
    print(f"{'Flexural Compression:':<{spacing}} {flexuralCompression(shapes)}")
    print(f"{'Maximum Shear:':<{spacing}} {centroidShear(shapes)}")
    print(f"{'Glue Shear:':<{spacing}} {glueShear(shapes, glue)}")
    c1, c2, c3 = splitCases(shapes)
    c1t, c1f = case1(c1, I, ybar)
    print(f"{'Case 1 Failure Stress:':<{spacing}} {min(c1t)}")
    c2t, c2f = case2(c2, I, ybar)
    print(f"{'Case 2 Failure Stress:':<{spacing}} {min(c2t)}")
    c3t, c3f = case3(c3, I, ybar)
    print(f"{'Case 3 Failure Stress:':<{spacing}} {min(c3t)}")
    print(f"{'Case 4 Failure Stress:':<{spacing}} {case4(shapes)}")
    print("-" * 80)

    print("\n" + "-" * 80)
    print(f"{'FOS Tension:':<{spacing}} {30 / flexuralTension(shapes)}")
    print(f"{'FOS Compression:':<{spacing}} {6 / flexuralCompression(shapes)}")
    print(f"{'FOS Glue:':<{spacing}} {2/glueShear(shapes, glue)}")
    print(f"{'FOS Shear:':<{spacing}} {min(4, case4(shapes))/centroidShear(shapes)}")
    print(f"{'FOS Case 1:':<{spacing}} {min(c1f)}")
    print(f"{'FOS Case 2:':<{spacing}} {min(c2f)}")
    print(f"{'FOS Case 3:':<{spacing}} {min(c3f)}")

    print("-" * 80)
    print("\n")
    print("\n")


    draw_cross_section(shapes)

#Outputs fos' such that they can be copy pasted into a csv or excel file and keep their tabular format
def outputCopyPaste(shapes):
    shapes = convert(shapes)
    c1, c2, c3 = splitCases(shapes)
    ybar = centroid(shapes)
    I = secondmomentarea(shapes, ybar)
    c1t, c1f = case1(c1, I, ybar)
    c2t, c2f = case2(c2, I, ybar)
    c3t, c3f = case3(c3, I, ybar)
    
    
    print(f"Height,{height(shapes)}")
    print(f"Centroid,{centroid(shapes)}")
    print(f"I,{I}")
    print(f"Diaphragms, {diaphragms}")

    print(f"Tension,{30 / flexuralTension(shapes)}")
    print(f"Compression,{6 / flexuralCompression(shapes)}")
    print(f"Glue,{2 / glueShear(shapes, glue)}")
    print(f"Shear,{min(4, case4(shapes)) / centroidShear(shapes)}")
    print(f"Case 1,{min(c1f)}")
    print(f"Case 2,{min(c2f)}")
    print(f"Case 3,{min(c3f)}")




    



if __name__ == "__main__":
    #outputCopyPaste(shapes)
    bundledOutput(shapes)

    # shapes = convert(shapes)
    # draw_cross_section(shapes)
    # print(f"Centroid: {centroid(shapes)} from y = 0")
    # print(f"I: {secondmomentarea(shapes, centroid(shapes))}")
    # print(glueShear(shapes, glue))
    # print(flexuralCompression(shapes))
    # print(flexuralTension(shapes))
    #print(centroidShear(shapes))
    # c1, c2, c3 = splitCases(shapes)
    # print(case1(c1))
    # print(case2(c2))
    # print(case3(c3))
    #print(case4(shapes))