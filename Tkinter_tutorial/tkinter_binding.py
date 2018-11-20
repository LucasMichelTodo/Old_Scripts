#!/usr/bin/env python

import Tkinter as tkinter

window = tkinter.Tk()
window.title("GUI")

# creating a function with an arguments 'event'
def say_hi(event): # you can rename 'event' to anything you want
    tkinter.Label(window, text = "Hi").pack()

btn = tkinter.Button(window, text = "Click Me!")
btn.bind("<Button-1>", say_hi) # 'bind' takes 2 parameters 1st is 'event' 2nd is 'function'
btn.pack()

window.mainloop()
