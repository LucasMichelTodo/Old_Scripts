#!/usr/bin/env python

import Tkinter as tkinter

class GeeksBro:

    def __init__(self, window):

        self.text_btn = tkinter.Button(window, text = "Click Me!", command = self.say_hi) # create a button to call a function called 'say_hi'
        self.text_btn.pack()

        self.close_btn = tkinter.Button(window, text = "Close", command = window.quit) # closing the 'window' when you click the button
        self.close_btn.pack()

    def say_hi(self):
        tkinter.Label(window, text = "Hi").pack()

window = tkinter.Tk()
window.title("GUI")

geeks_bro = GeeksBro(window)

window.mainloop()
