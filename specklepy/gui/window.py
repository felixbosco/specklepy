import sys
from PyQt5 import QtWidgets, QtGui, QtCore

from specklepy.logging import logger


def start():
    logger.info("Initializing application...")
    app = QtWidgets.QApplication(sys.argv)

    logger.info("Initializing GUI...")
    win = Window()
    win.show()

    sys.exit(app.exec_())


class Window(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()

        # Settings
        self.setGeometry(50, 50, 1000, 600)
        self.setWindowTitle('Specklepy')
        # self.setWindowIcon(QtGui.QIcon('source/img/telescope.png'))

        # Menu bar
        main_menu = self.menuBar()
        self.create_dropdown(main_menu, '&File', [{'text': 'New window', 'function': self.__init__,
                                                   'shortcut': 'Ctrl+Shift+N'},
                                                  {'text': 'Open', 'function': self.do_nothing, 'shortcut': 'Ctrl+O'},
                                                  {'text': 'Settings', 'function': self.settings},
                                                  {'text': 'Quit', 'function': self.close_application,
                                                  'shortcut': 'Ctrl+Q'}])
        self.create_dropdown(main_menu, '&Edit', [{'text': 'Undo', 'function': self.do_nothing,
                                                   'status_tip': 'Do nothing'},
                                                  {'text': 'Redo', 'function': self.do_nothing,
                                                  'status_tip': 'Do nothing'}])
        self.create_dropdown(main_menu, '&View', [{'text': 'Increase font size', 'function': self.do_nothing,
                                                   'shortcut': 'Ctrl++'},
                                                  {'text': 'Decrease font size', 'function': self.do_nothing,
                                                  'shortcut': 'Ctrl+-'}])
        self.create_dropdown(main_menu, '&Help', [{'text': 'Call help', 'function': self.do_nothing,
                                                   'shortcut': 'Ctrl+H'}])

        # Background
        self.setAutoFillBackground(True)
        palette = self.palette()
        palette.setColor(self.backgroundRole(), QtCore.Qt.white)
        self.setPalette(palette)

        # View
        self.home()

    # Views
    def home(self):
        btn = self.create_button('Quit', self.close_application, move=[0, 500])
        btn = self.create_button('Settings', self.settings, move=[200, 200])
        self.show()

    def settings(self):
        btn = self.create_button('Quit', self.close_application, move=[0, 500])
        toolBar = self.addToolBar('Tools')
        toolBar.addAction(self.create_action('Quit', self.close_application, icon='source/img/squares_100x100.png'))
        self.show()

    # Methods
    def create_menu(self, dropdowns):
        menu = self.menuBar()
        for dropdown in dropdowns:
            self.create_dropdown(menu, **dropdown)
        return menu

    def create_dropdown(self, menu, text, actions):
        dropdown = menu.addMenu(text)
        for action in actions:
            dropdown.addAction(self.create_action(**action))
        return dropdown

    def close_application(self):
        choice = QtWidgets.QMessageBox.question(self, 'Close', 'Do you really want to close the application?',
                                                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            logger.info("Exiting application...")
            sys.exit()

    def create_action(self, text, function, **kwargs):
        if 'icon' in kwargs:
            action = QtWidgets.QAction(QtGui.QIcon(kwargs['icon']), text, self)
        else:
            action = QtWidgets.QAction(text, self)
        action.triggered.connect(function)
        if 'shortcut' in kwargs:
            action.setShortcut(kwargs['shortcut'])
        if 'status_tip' in kwargs:
            action.setStatusTip(kwargs['status_tip'])
        return action

    def create_button(self, text, function, **kwargs):
        btn = QtWidgets.QPushButton(text, self)
        btn.clicked.connect(function)
        if 'resize' in kwargs:
            btn.resize(kwargs['resize'])
        else:
            btn.resize(btn.sizeHint())
        if 'move' in kwargs:
            btn.move(*kwargs['move'])
        return btn

    def do_nothing(self):
        pass
