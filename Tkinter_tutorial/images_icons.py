#!/usr/bin/env python

import Tkinter as tkinter

window = tkinter.Tk()
window.title("GUI")

# taking image from the directory and storing the source in a variable
icon = tkinter.PhotoImage(file = "images/haha.png")
# displaying the picture using a 'Label' by passing the 'picture' variriable to 'image' parameter
label = tkinter.Label(window, image = icon)
label.pack()

window.mainloop()
