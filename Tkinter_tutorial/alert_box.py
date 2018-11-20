#!/usr/bin/env python

import Tkinter as tkinter


window = tkinter.Tk()
window.title("GUI")

# creating a simple alert box
tkinter.messagebox.showinfo("Alert Message", "This is just a alert message!")
# creating a question to get the response from the user [Yes or No Question]
response = tkinter.messagebox.askquestion("Simple Question", "Do you love Python?")
# If user clicks 'Yes' then it returns 1 else it returns 0
if response == 1:
    tkinter.Label(window, text = "You love Python!").pack()
else:
    tkinter.Label(window, text = "You don't love Python!").pack()

window.mainloop()
