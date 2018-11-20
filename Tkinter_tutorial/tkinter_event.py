c

## Click and say Hi!

window = tkinter.Tk()
window.title("GUI")

# creating a function called say_hi()
def say_hi():
    tkinter.Label(window, text = "Hi").pack()

tkinter.Button(window, text = "Click Me!", command = say_hi).pack() # 'command' is executed when you click the button
                                                                    # in this above case we're calling the function 'say_hi'.
window.mainloop()
