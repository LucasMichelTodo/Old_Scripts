#!/usr/bin/env python

import Tkinter as tkinter

window = tkinter.Tk()
window.title("GUI")

# creating the 'Canvas' area of width and height 500px
canvas = tkinter.Canvas(window, width = 500, height = 500)
canvas.pack()

# 'create_line' is used to create a line. Parameters:- (starting x-point, starting y-point, ending x-point, ending y-point)
line1 = canvas.create_line(25, 25, 250, 150)
# parameter:- (fill = color_name)
line2 = canvas.create_line(25, 250, 250, 150, fill = "red")

# 'create_rectangle' is used to create rectangle. Parameters:- (starting x-point, starting y-point, width, height, fill)
# starting point the coordinates of top-left point of rectangle
rect = canvas.create_rectangle(500, 25, 175, 75, fill = "green")

# you 'delete' shapes using delete method passing the name of the variable as parameter.
#canvas.delete(line1)
# you 'delete' all the shapes by passing 'ALL' as parameter to the 'delete' method
# canvas.delete(tkinter.ALL)

window.mainloop()
