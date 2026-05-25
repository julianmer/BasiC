"""
main.py
Basis Set Converter — entry point

Launches the app window and manages navigation between screens.
Each screen is its own class in its own file.

Shared state is a plain dict passed to every screen.
Screens read from it and write to it as needed.

State keys (added progressively as the user moves through screens):
    'files'            : list of file paths uploaded by the user
    'detected_format'  : string, e.g. 'LCModel .BASIS'
    'detection_status' : 'full' | 'partial' | 'ambiguous' | 'failed'
    'core'             : dict — the core struct (fids, sw, sf, n, names, ...)
    'missing_fields'   : list of field names that could not be found or derived
    'output_tool'      : string, e.g. 'Osprey'
    'output_format'    : string, e.g. 'mat'
    'output_path'      : string, path where output was saved
    'params'           : dict — user-supplied or defaulted parameter values
"""

import tkinter as tk

from screen1_upload  import Screen1Upload
from screen2_output  import Screen2Output
from screen3_params  import Screen3Params
from screen4_convert import Screen4Convert


class App(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Basis Set Converter")
        self.resizable(True, True)
        self.minsize(700, 500)

        # Shared state dict — passed to every screen
        self.state = {}

        # Build all screens upfront, stacked in the same grid cell
        self.screens = {}
        self._build_screens()

        # Start on screen 1
        self.show_screen("screen1")

    def _build_screens(self):
        """Instantiate all screen frames and place them in the window."""
        screen_classes = {
            "screen1": Screen1Upload,
            "screen2": Screen2Output,
            "screen3": Screen3Params,
            "screen4": Screen4Convert,
        }

        for name, cls in screen_classes.items():
            frame = cls(self, self.state, self)
            frame.grid(row=0, column=0, sticky="nsew")
            self.screens[name] = frame

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

    def show_screen(self, name):
        """
        Show the named screen, hide all others.
        Called by screens to navigate forward or backward.

        Usage from any screen:
            self.app.show_screen("screen2")
        """
        for screen_name, frame in self.screens.items():
            if screen_name == name:
                frame.tkraise()
                frame.on_show()


if __name__ == "__main__":
    app = App()
    app.mainloop()
