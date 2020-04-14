# -*- coding: utf-8 -*-
"""
Created on Tue May 24 13:09:08 2016

@author: Josh Bradshaw - joshbradshaw11@gmail.com
"""
from __future__ import absolute_import, print_function, division

# std lib imports
import os
from functools import wraps
import json
import math

# anaconda module imports
import qtpy
from qtpy import QtCore, QtWidgets
import qtawesome
if qtpy.API == 'pyqt5':
    from matplotlib.backends.backend_qt5agg import \
        FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
else:
    from matplotlib.backends.backend_qt4agg import \
        FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
# other python files that should be in the same dir
from . import ROI
from . import blood_tools
from . import fitting


def callback(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        try:
            return f(*args, **kwds)
        except Exception as e:
            print(e)
    return wrapper


class ROISelectPlot(QtWidgets.QWidget):
    """
    The plot canvas for the image, in the left and center panes of the GUI. 
    This is where the image is displayed, and the ROI is selected.
    """

    def __init__(self, parent=None):
        super(ROISelectPlot, self).__init__(parent)
        # initialize the plot area
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)
        self.im = None
        self.mpl_im = None

    def make_image(self, im, colormap, vmin=5, vmax=95):
        self.im = im
        ax = self.axes
        self.mpl_im = ax.imshow(
            im, vmin=np.percentile(im, vmin), vmax=np.percentile(im, vmax),
            cmap=colormap, origin='upper',
            extent=[0, self.im.shape[1], self.im.shape[0], 0])
        ax.set_xlim(0, self.im.shape[1])
        ax.set_ylim(self.im.shape[0],0)
        self.toolbar.update()
        self.toolbar.push_current()
        self.figure.canvas.draw()


class T2CurvePlot(QtWidgets.QWidget):
    """
    The plot that displays the datapoints, the fitted monoexponential T2 curve, 
    and gives the value for T2 based on this fit.
    """

    def __init__(self, parent=None):
        super(T2CurvePlot, self).__init__(parent)
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111, autoscale_on=True)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)


class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        # draw the interface
        self.init_image_data()
        self.init_gui()

        # init controls
        self.vmin = 5
        self.vmax = 95
        self.controls_enabled(False)

    def init_image_data(self):
        # initialize these to empty lists so that the slice select display will
        # work before loading images
        self.images, self.image_attributes = [], []
        # used for keeping track of what slice is displayed
        self.image_index = 0  # index of the image that is currently displayed
        self.image_filename = None  # filename of the image that is currently displayed
        # dictionary that stores image ROIs, the key is the image filename
        # I opted to do it this way so that the ROIs can be reloaded
        self.image_ROIs = {}
        self.image_filename_list = []  # list of image filenames for the current series
        self.grey_activeROI = None  # ROI object displayed on the greyscale plot
        self.color_activeROI = None  # ROI object displayed on the colour plot
        self.roi_path = None  # filename where ROI data is saved
        self.directory = ""  # DICOM directory, selected by the user
        # matplotlib object generated when the ROI is drawn on the image
        self.grey_roi_patch = None
        # matplotlib object generated when the ROI is drawn on the image
        self.color_roi_patch = None
        self.included_slices = []

    def init_gui(self):
        QtWidgets.QWidget.__init__(self, parent=None)
        # button for opening dicom directory
        self.button_load = QtWidgets.QPushButton(
            qtawesome.icon('fa.folder-open-o'), '')
        self.button_run = QtWidgets.QPushButton('Fit Data')
        self.button_save_ROI = QtWidgets.QPushButton('Save ROIs')
        self.button_load_ROI = QtWidgets.QPushButton('Load ROIs')
        self.button_draw_roi = QtWidgets.QPushButton('Draw ROI')
        self.button_exclude_slice = QtWidgets.QPushButton('Exclude Slice')
        # first/last prev/next buttons for scrolling through the dicoms
        self.button_image_fwd = QtWidgets.QPushButton(
            qtawesome.icon('fa.chevron-right'), '')
        self.button_image_fwd.setObjectName('slice_fwd')
        self.button_image_bwd = QtWidgets.QPushButton(
            qtawesome.icon('fa.chevron-left'), '')
        self.button_image_bwd.setObjectName('slice_bwd')

        self.button_image_first = QtWidgets.QPushButton(
            qtawesome.icon('fa.step-backward'), '')
        self.button_image_first.setObjectName('slice_first')
        self.button_image_last = QtWidgets.QPushButton(
            qtawesome.icon('fa.step-forward'), '')
        self.button_image_last.setObjectName('slice_last')
        self.slice_label = QtWidgets.QLabel('(00/00)')
        # set whether the ROI should apply to a single slice or all of the
        # slices
        self.combo_roi_scope = QtWidgets.QComboBox()
        self.combo_roi_scope.addItems(['This Slice', 'All Slices'])
        self.combo_roi_scope.setCurrentIndex(0)
        # choose ROI shape
        self.combo_roi_style = QtWidgets.QComboBox()
        self.combo_roi_style.addItems(['Polygon', 'Circle', 'Ellipse'])
        self.combo_roi_style.setCurrentIndex(0)
        # choose colormap
        self.combo_colormap = QtWidgets.QComboBox()
        self.combo_colormap.addItems(['gray', 'jet', 'hot', 'viridis'])
        self.combo_colormap.setCurrentIndex(0)
        # fit either a basic T1, or basic T2 fit
        self.combo_relax_label = QtWidgets.QLabel('Fit Type')
        self.combo_relax = QtWidgets.QComboBox()
        self.combo_relax.addItems(['T1', 'T2'])
        self.combo_relax.setCurrentIndex(0)
        
        self.combo_type_label = QtWidgets.QLabel('Fit parameters')
        self.combo_type = QtWidgets.QComboBox()
        self.combo_type.addItems(['2', '3'])
        self.combo_type.setCurrentIndex(0)
        

        self.roi_area_label = QtWidgets.QLabel(
            'ROI Area: 0.00 (pixels) / 0.00 (mm^2)')
        # doesn't actually do anything
       # self.uncertainty_checkbox = QtWidgets.QCheckBox('Fit with uncertainty', self)

        self.plot_im = ROISelectPlot(self)
        self.plot_graph = T2CurvePlot(self)

        layout_top = QtWidgets.QHBoxLayout()
        layout_top.addSpacing(10)
        layout_top.addWidget(self.button_load)

        # layout_top.addStretch()
        layout_top.addWidget(QtWidgets.QLabel('Change Slice:'))

        layout_top.addWidget(self.button_image_first)
        layout_top.addWidget(self.button_image_bwd)
        layout_top.addWidget(self.button_image_fwd)
        layout_top.addWidget(self.button_image_last)
        layout_top.addWidget(self.slice_label)

        layout_top.addStretch()

        layout_subTop = QtWidgets.QHBoxLayout()
        layout_subTop.addWidget(QtWidgets.QLabel('ROI Style:'))
        layout_subTop.addWidget(self.combo_roi_style)
        layout_subTop.addWidget(QtWidgets.QLabel('Apply ROI to:'))
        layout_subTop.addWidget(self.combo_roi_scope)
        layout_subTop.addWidget(self.button_draw_roi)
        layout_subTop.addWidget(self.button_exclude_slice)
        layout_subTop.addWidget(self.button_save_ROI)
        layout_subTop.addWidget(self.button_load_ROI)
        layout_subTop.addStretch()

        layout_sub_subTop = QtWidgets.QHBoxLayout()
        layout_sub_subTop.addStretch()
        layout_sub_subTop.addWidget(self.combo_type_label)
        layout_sub_subTop.addWidget(self.combo_type)
        layout_sub_subTop.addWidget(self.combo_relax_label)
        layout_sub_subTop.addWidget(self.combo_relax)
        # layout_top.addWidget(self.uncertainty_checkbox)
        layout_sub_subTop.addWidget(self.button_run)

        layout_mid = QtWidgets.QHBoxLayout()
        layout_mid.addWidget(self.plot_im)
        layout_mid.addWidget(self.plot_graph)

        self.vmin_window_slider = QtWidgets.QSlider(
            orientation=QtCore.Qt.Horizontal, minimum=0, maximum=100)
        self.vmin_window_slider.setValue(5)
        self.vmax_window_slider = QtWidgets.QSlider(
            orientation=QtCore.Qt.Horizontal, minimum=0, maximum=100)
        self.vmax_window_slider.setValue(95)

        layout_ROI_calc = QtWidgets.QHBoxLayout()
        layout_ROI_calc.addWidget(QtWidgets.QLabel('Colormap:'))
        layout_ROI_calc.addWidget(self.combo_colormap)
        layout_ROI_calc.addStretch()
        layout_ROI_calc.addWidget(self.roi_area_label)
        layout_slider1 = QtWidgets.QHBoxLayout()
        layout_slider1.addWidget(QtWidgets.QLabel('Window Min:'))
        layout_slider1.addSpacing(3)
        layout_slider1.addWidget(self.vmin_window_slider)

        layout_slider2 = QtWidgets.QHBoxLayout()
        layout_slider2.addWidget(QtWidgets.QLabel('Window Max:'))
        layout_slider2.addWidget(self.vmax_window_slider)

        layout_main = QtWidgets.QVBoxLayout()
        layout_main.addLayout(layout_top)
        layout_main.addLayout(layout_subTop)
        layout_main.addLayout(layout_sub_subTop)
        layout_main.addLayout(layout_mid)
        layout_main.addLayout(layout_ROI_calc)
        layout_main.addLayout(layout_slider1)
        layout_main.addLayout(layout_slider2)
        self.setLayout(layout_main)

        self.button_load.pressed.connect(self.choose_dir)
        self.button_run.pressed.connect(self.process_data)
        self.button_draw_roi.pressed.connect(self.start_roi)
        self.button_exclude_slice.pressed.connect(self.exclude_current_slice)
        self.button_image_first.pressed.connect(self.change_image)
        self.button_image_last.pressed.connect(self.change_image)
        self.button_image_fwd.pressed.connect(self.change_image)
        self.button_image_bwd.pressed.connect(self.change_image)
        self.button_save_ROI.pressed.connect(self.save_ROI)
        self.button_load_ROI.pressed.connect(self.load_prev_analysis)
        self.vmin_window_slider.valueChanged.connect(self.set_image_window)
        self.vmax_window_slider.valueChanged.connect(self.set_image_window)
        self.combo_colormap.currentIndexChanged.connect(self.colormap_changed)

    @callback
    def colormap_changed(self, *e):
        colormap = self.combo_colormap.currentText()
        self.plot_im.make_image(
            self.images[self.image_index], colormap, self.vmin, self.vmax)

    @callback
    def exclude_current_slice(self, *e):
        """allows to exclude a motion corrupted slice"""
        if self.included_slices[self.image_index]:
            self.included_slices[self.image_index] = False
            self.clear_roi()
            self.roi_controls_enable(False)
        else:
            self.included_slices[self.image_index] = True
            self.roi_controls_enable(True)
            self.load_roi()

    @callback
    def save_ROI(self, *e):
        to_save = {}
        image_ROIs_json = self.create_image_ROIs_json()
        # to_save['ROIs'] = self.image_ROIs
        to_save['ROIs'] = image_ROIs_json
        to_save['included_slices'] = self.included_slices
        if self.directory:
            out = QtWidgets.QFileDialog.getSaveFileName(
                directory=self.directory, caption='ROI filename')
        else:
            out = QtWidgets.QFileDialog.getSaveFileName(caption='ROI filename')
        if len(out) == 2:
            out = out[0]
        if not out:
            return
        if not out.endswith('.json'):
            out = out + '.json'
        self.roi_path = out
        with open(self.roi_path, 'w') as f:
            json.dump(to_save, f)

    def roi_controls_enable(self, enable=True):
        """disable the ROI controls when the slice is excluded, to provide
        a salient visual cue that this slice isn't being taken into account"""
        if enable:
            self.button_exclude_slice.setText('Exclude Slice')
        else:
            self.button_exclude_slice.setText('Include Slice')

        self.combo_roi_style.setEnabled(enable)
        self.combo_relax.setEnabled(enable)
        self.combo_type.setEnabled(enable)
        self.combo_roi_scope.setEnabled(enable)
        self.button_draw_roi.setEnabled(enable)

    def controls_enabled(self, enable=True):
        """disable or enable all controls, to keep the user from pressing buttons
        that have no function when datasets aren't yet loaded"""
        self.button_run.setEnabled(enable)
        self.button_draw_roi.setEnabled(enable)
        self.button_image_first.setEnabled(enable)
        self.button_image_last.setEnabled(enable)
        self.button_image_fwd.setEnabled(enable)
        self.button_image_bwd.setEnabled(enable)
        self.button_exclude_slice.setEnabled(enable)
        self.combo_roi_style.setEnabled(enable)
        self.combo_colormap.setEnabled(enable)
        self.combo_relax.setEnabled(enable)
        self.combo_type.setEnabled(enable)
        self.combo_roi_scope.setEnabled(enable)
        self.button_save_ROI.setEnabled(enable)
        self.button_load_ROI.setEnabled(enable)

    @callback
    def set_image_window(self, *e):
        self.vmin = self.vmin_window_slider.value()  # image window minimum
        self.vmax = self.vmax_window_slider.value()  # image window maximum
        
        self.plot_im.make_image(
            self.images[self.image_index], self.combo_colormap.currentText(), self.vmin, self.vmax)

    def clear_roi(self):
        # remove the ROI from the screen, but do not delete it until it is
        # overwritten by another ROI
        # remove the ROI patches created when loading if necessary
        if self.grey_roi_patch is not None:
            self.grey_roi_patch.remove()
            self.grey_roi_patch = None
        if self.color_roi_patch is not None:
            self.color_roi_patch.remove()
            self.color_roi_patch = None

        if self.grey_activeROI is not None:  # check if there is an ROI
            self.grey_activeROI.remove()
            self.grey_activeROI = None

        grey_figure = self.plot_im.figure
        grey_figure.canvas.draw()

    def change_image(self):
        # serialize existing ROIs to file, this is quick+dirty b/c I haven't
        # figured out how to detect the ROI complete event in this class
        # so I just save them when the user changes images
        sender_btn = self.sender().objectName()  # <<, <, >, >> are the possibilities
        num_images = len(self.images)

        # prevent the buttons from raising div by zero exceptions when no
        # images loaded
        if num_images > 0:
            if sender_btn == 'slice_first':
                if self.t2:
                    self.image_index = np.argsort(self.prep_times)[0]
                else:
                    self.image_index = 0
            elif sender_btn == 'slice_bwd':
                if self.t2:
                    ind = np.where(np.argsort(self.prep_times)
                                   == self.image_index)[0][0] - 1
                    self.image_index = np.argsort(self.prep_times)[
                        ind % num_images]
                else:
                    ind = self.image_index - 1
                    self.image_index = ind % num_images
            elif sender_btn == 'slice_fwd':
                if self.t2:
                    ind = np.where(np.argsort(self.prep_times)
                                   == self.image_index)[0][0] + 1
                    self.image_index = np.argsort(self.prep_times)[
                        ind % num_images]
                else:
                    ind = self.image_index + 1
                    self.image_index = ind % num_images
            else:
                if self.t2:
                    ind = num_images - 1
                    self.image_index = np.argsort(self.prep_times)[ind]
                else:
                    self.image_index = num_images - 1

            self.image_filename = self.image_filename_list[self.image_index]
            # display the slice selection label, with zero padding to keep the
            # toolbar from shifting around
            num, demon = str(self.image_index + 1).rjust(
                2, '0'), str(num_images).rjust(2, '0')
            # display previous ROI if it exists
            self.clear_roi()
            self.plot_im.mpl_im.set_data(self.images[self.image_index])
            axes = self.plot_im.axes
            if self.t2:
                title_str = 'TE=%.0f ms' % self.prep_times[self.image_index]
            else:
                inversion_times = [att['InversionTime']
                                   for att in self.image_attributes]
                title_str = 'TI=%.0f ms' % inversion_times[self.image_index]
            axes.set_title(title_str)

            self.plot_im.figure.canvas.draw()
            self.load_roi()
            self.slice_label.setText("{}/{}".format(num, demon))

    def create_image_ROIs_json(self):
        """create a version of the self.image_ROIs dictionary where the keys are json strings instead of ROI objects"""
        image_ROIs_json = {}
        for img_fn in self.image_filename_list:
                image_ROIs_json[img_fn]=ROI.write_json_str(self.image_ROIs[img_fn])
        return image_ROIs_json
    
    def create_image_ROIs_from_json(self, ROIs_json):
        """Retrieve original version of image_ROIs from dictionary with json strings"""
        image_ROIs_from_json = {}
        for img_fn in self.image_filename_list:
                image_ROIs_from_json[img_fn] = ROI.load_from_json(
                    ROIs_json[img_fn])
        return image_ROIs_from_json 

    @callback
    def load_prev_analysis(self):
        """Load the previous ROIs and slice inclusion list. I opted to prompt the user
        so that they understand where the ROIs come from as opposed to having them
        magically appear
        """
        self.clear_roi()
        if self.directory:
            infilename = QtWidgets.QFileDialog.getOpenFileName(
                directory=self.directory, caption='ROI filename')
        else:
            infilename = QtWidgets.QFileDialog.getOpenFileName(
                caption='ROI filename')
        if len(infilename) == 2:
            infilename = infilename[0]
        if not infilename:
            return
        with open(infilename, 'r') as f:
            to_load = json.load(f)
        self.image_ROIs = self.create_image_ROIs_from_json(to_load['ROIs'])
        self.included_slices = to_load['included_slices']
        self.load_roi()

    def calc_ROI_area(self):
        """using pixel spacing tag in the DICOM header, calculate the area
        represented by the pixels in the ROI

        NOTE: THIS MAY NOT WORK PROPERLY, TODO: TEST mm^2 area
        """
        indicies = self.image_ROIs[self.image_filename].get_indices()
        pixel_count = len(indicies[0])
        row_spacing, col_spacing = self.image_pixel_spacing
        area_mm2 = pixel_count * row_spacing * col_spacing

        roi_area_str = 'ROI Area: {0} (pixels) / {1:.2f} (mm^2)'.format(
            pixel_count, area_mm2)

        self.roi_area_label.setText(roi_area_str)

    # todo: rework these into one simpler method (its tricky because they're
    # callbacks!)
    def grey_roi_complete_callback(self):
        """method called when user finishes drawing an ROI on the greyscale plot
        mirrors the ROI over the colour plot and saves the ROI"""
        roi_scope = self.combo_roi_scope.currentText()

        if roi_scope == "All Slices":
            for img_fn in self.image_filename_list:
                self.image_ROIs[img_fn] = self.grey_activeROI
        else:
            self.image_ROIs[self.image_filename] = self.grey_activeROI

        self.clear_roi()
        self.load_roi()

    @callback
    def choose_dir(self, *event):
        """opens a directory choose dialog box, allows the user to select their
        dicom series of interest and loads that series."""
        self.plot_graph.axes.clear()
        self.plot_graph.figure.canvas.draw()

        if self.directory:
            out = QtWidgets.QFileDialog.getExistingDirectory(
                directory=os.path.split(self.directory)[0], caption='MRI Data Directory')
        else:
            out = QtWidgets.QFileDialog.getExistingDirectory(
                caption='MRI Data Directory')

        if out:
            self.image_index = 0
            self.t2 = 0  # is it a T2 series?
            self.image_filename = None
            self.image_ROIs = {}
            self.image_filename_list = []
            self.grey_activeROI = None
            self.color_activeROI = None
            self.roi_path = None
            self.images, self.image_attributes = [], []
            self.plot_graph.axes.clear()
            self.clear_roi()
            self.directory = out
            self.roi_path = os.path.join(self.directory, '.ROIs')
            self.images, self.image_attributes, self.dicom_list = blood_tools.read_dicoms(
                out, ['InversionTime', 'PixelSpacing'])

            VB17_prep_times = blood_tools.get_T2_prep_times_VB17(
                self.dicom_list)
            VE11_prep_times = blood_tools.get_T2_prep_times_VE11(
                self.dicom_list)
            VD13_prep_times = blood_tools.get_T2_prep_times_VD13A(
                self.dicom_list)

            if VE11_prep_times:
                self.prep_times = VE11_prep_times
            if VB17_prep_times:
                self.prep_times = VB17_prep_times
            if VD13_prep_times:
                self.prep_times = VD13_prep_times
                
            #for some reason, zero prep time image gets no prep time, assign it TE=2000    
            if len(self.prep_times) < len(self.dicom_list):
                 self.prep_times=np.concatenate([np.array([2000]),self.prep_times])       

            if not self.images:
                error = QtWidgets.QErrorMessage()
                error.showMessage(
                    'The selected directory does not contain a DICOM series '
                    'which this widget is capable of loading')
                error.exec_()
                return

            if len(self.prep_times) == len(self.images):
                self.t2 = 1
                self.image_index = np.argsort(self.prep_times)[0]
                self.combo_relax.setCurrentIndex(1)
            self.image_pixel_spacing = self.image_attributes[0]['PixelSpacing']
            # initiate included slices to be all True
            self.included_slices = [True for _ in range(len(self.images))]

            for attributes in self.image_attributes:
                self.image_filename_list.append(attributes['filename'])

            self.image_filename = self.image_filename_list[self.image_index]
            colormap = self.combo_colormap.currentText()
            self.plot_im.make_image(
                self.images[self.image_index], colormap, self.vmin, self.vmax)
            axes = self.plot_im.axes
            if self.t2:
                title_str = 'TE=%.0f ms' % self.prep_times[self.image_index]
            else:
                inversion_times = [att['InversionTime']
                                   for att in self.image_attributes]
                title_str = 'TI=%.0f ms' % inversion_times[self.image_index]
            axes.set_title(title_str)
            self.controls_enabled(True)
            num, demon = '01', str(len(self.images)).rjust(2, '0')
            self.slice_label.setText("{}/{}".format(num, demon))

        else:  # user hit the cancel or x button to leave the dialog
            pass

    def load_roi(self):
        """if the active image has a previously drawn ROI, this method reloads its"""

        if not self.included_slices[self.image_index]:
            self.roi_controls_enable(False)
            return
        else:
            self.roi_controls_enable(True)

        if self.image_filename in self.image_ROIs:
            self.grey_activeROI = self.image_ROIs[self.image_filename]
            self.color_activeROI = self.image_ROIs[self.image_filename]
            self.grey_roi_patch = self.grey_activeROI.draw(
                self.plot_im.axes, self.plot_im.figure, 'red')
            self.calc_ROI_area()

    @callback
    def start_roi(self):
        """create a new ROI for the image"""
        if not len(self.images) > 0:
            error = QtWidgets.QErrorMessage()
            error.showMessage(
                'You must a load a series of images before drawing the ROI')
            error.exec_()
            return

        roi_style = self.combo_roi_style.currentText().lower()
        # style names are lowercase in ROI.py

        self.clear_roi()
        # create an ROI object for both images, keep the one that calls the
        # complete callback first
        self.grey_activeROI = ROI.new_ROI(self.plot_im.mpl_im,
                                          self.plot_im.axes,
                                          self.plot_im.figure,
                                          roi_style, 'red',
                                          self.grey_roi_complete_callback)

    @callback
    def process_data(self, *event):
        """Gets the prep times and populates the T1/T2 plot"""
        # check that user has drawn all of the required ROIs
        for ii, fn in enumerate(self.image_filename_list):
            if self.included_slices[ii] and not fn in self.image_ROIs:
                error = QtWidgets.QErrorMessage()
                error.showMessage(
                    "You must draw an ROI on every included slice before you can fit the data. Use the 'All Slices' option if the ROIs are conincident accross the slices")
                error.exec_()
                return

        roi_list = []
        for image_filename in self.image_filename_list:
            roi_list.append(self.image_ROIs[image_filename])

        if not len(self.images) > 0:
            error = QtWidgets.QErrorMessage()
            error.showMessage(
                'You must load a series of dicom images before fitting the data')
            error.exec_()
            return

        # todo add error message if images or ROI not loaded
        relaxation_type = self.combo_relax.currentText()
        axes = self.plot_graph.axes
        axes.clear()
        
        fit_type=self.combo_type.currentText()

        if relaxation_type == 'T1':
            self.t2=0
            ti, signal, stddev = get_T1_decay_signal(
                self.image_attributes, self.images, roi_list, self.included_slices)
            
            if not len(ti) or ti[0] == 0:
                error = QtWidgets.QErrorMessage()
                error.showMessage(
                    'Failed to find T1 recovery times for this dataset, ensure that this is a T1 series')
                error.exec_()
                return

            axes.errorbar(ti, signal, stddev, fmt='o')
            axes.set_xlabel('inversion time (ms)')
            axes.set_ylabel('signal intensity')


#            inversion_recovery = fitting.model(
#                'abs(M0*(1-2*aa*exp(-x/T1)))', {'M0': signal.max(), 'aa': 1, 'T1': 1000})
            
            inversion_recovery = fitting.model(
                    'abs(A-B*exp(-x/T1))',{'A':signal.max(),'B':0.5*signal.max(),'T1':1000})

            # Initial T1 guess based on null time
            ti = np.array(ti)  # necessary for the next line to work
            if 20000 not in ti:
                T1_guess = -ti[np.where(signal == signal.min())] / np.log(0.5)
                inversion_recovery['T1'] = T1_guess[0]        
                inversion_recovery['B']=2*signal.max()        

            bootstrap_signals = bootstrap_decay_signal(
                self.image_attributes, self.images, roi_list, self.included_slices, iterations=100)
            inversion_recovery.fit(ti, signal, stddev)
            par_uncertainties = fit_bootstrap(
                ti, bootstrap_signals, inversion_recovery)

            fix_x_points = np.arange(0, 1.3 * np.max(ti), 1)
            axes.plot(fix_x_points, inversion_recovery(fix_x_points), 'b')
#            T1_corr = inversion_recovery['T1'].value * \
#                (2 * inversion_recovery['aa'].value - 1)
            
            T1_corr = inversion_recovery['T1'].value*(inversion_recovery['B'].value/inversion_recovery['A'].value-1)

            # Use error propagation to obtain uncertainty in LL corrected T1
#            temp = np.sqrt((par_uncertainties['aa'] / inversion_recovery['aa'].value)**2 + (
#                par_uncertainties['T1'] / inversion_recovery['T1'].value)**2)
#            T1_corr_unc = np.sqrt(temp**2 + par_uncertainties['T1']**2)
            
            T1_corr_unc = T1_corr*np.sqrt((par_uncertainties['A']/inversion_recovery['A'].value)**2 +
                           (par_uncertainties['B']/inversion_recovery['B'].value)**2 + 
                           (par_uncertainties['T1']/inversion_recovery['T1'].value)**2);                      

            resids = (inversion_recovery(ti) - np.array(signal)
                      )**2 / np.array(signal)**2
            RMSD = 100 * np.sqrt(np.mean(resids))

            str1 = r'$\rm{T_1=%.0f \pm %.0f\ ms}$' % (
                inversion_recovery['T1'].value, par_uncertainties['T1'])

            if 20000 in ti:
                str2=''
            else:
                str2 = r'$\rm{T_1\ (LL\ corr)=%.0f \pm %.0f\ ms}$' % (
                    T1_corr, T1_corr_unc)
            
            str3 = r'$\rm{RMSD=%.1f \%%}$' % (RMSD)

            #axes.text(0.5, 0.5, "T1 value, LL corrected : {}ms".format(round(T1_corr),0), transform=axes.transAxes)
            axes.text(0.4, 0.2, str1, transform=axes.transAxes, fontsize=16)
            axes.text(0.4, 0.1, str2, transform=axes.transAxes, fontsize=16)
            axes.set_title(str3, fontsize=16)

        elif relaxation_type == 'T2':

            te, signal, stddev = get_T2_decay_signal(
                self.dicom_list, self.images, roi_list, self.included_slices)
            # throw an error if no T2 prep times found in the dicom header
            if not len(te):
                error = QtWidgets.QErrorMessage()
                error.showMessage(
                    'Failed to find prep times for this dataset, ensure that this is a T2 series')
                error.exec_()
                return

            axes.errorbar(te, signal, stddev, fmt='o')
            axes.set_xlabel('TE time (ms)')
            axes.set_ylabel('signal intensity')

            start, stop = 0.5 * np.min(te), 1.25 * np.max(te)
            fit_x_points = np.linspace(start, stop, 1000)
            
            if fit_type == '2':
                spin_echo = fitting.model(
                    'M0*exp(-x/T2)', {'M0': np.max(signal), 'T2': 150})
                
            elif fit_type=='3':
                 spin_echo = fitting.model(
                    'M0*exp(-x/T2) + B', {'M0': np.max(signal), 'T2': 150, 'B':signal[0]})  

            #spin_echo.fit(te, signal, stddev)
            spin_echo.fit(te, signal, np.ones_like(signal))
            axes.plot(fit_x_points, spin_echo(fit_x_points), 'b')

            resids = (spin_echo(np.array(te)) -
                      np.array(signal))**2 / np.array(signal)**2
            RMSD = 100 * np.sqrt(np.mean(resids))

            bootstrap_signals = bootstrap_decay_signal(
                self.image_attributes, self.images, roi_list, self.included_slices, iterations=100)
            par_uncertainties = fit_bootstrap(te, bootstrap_signals, spin_echo)
            str1 = r'$\rm{T_2=%.0f \pm %.0f\ ms, RMSD=%.1f \%%}$' % (
                spin_echo['T2'].value, par_uncertainties['T2'], RMSD)
            #axes.text(0.3, 0.9, str1, transform=axes.transAxes, fontsize=16)
            axes.set_title(str1, fontsize=16)
        else:
            raise NotImplementedError(
                "This mapping type's fitting algorithm has not been implemented yet")
        axes.set_xlim(xmin=-5)
        self.plot_graph.figure.canvas.draw()


def get_T2_decay_signal(dicom_list, image_list, roi_list, included_slices, log_scale=False):
    """Gets T2 preps and signal"""
    # get T2 prep times
    # VE11 is the new Siemens software, VE17 is the old version
    VB17_prep_times = blood_tools.get_T2_prep_times_VB17(dicom_list)
    VE11_prep_times = blood_tools.get_T2_prep_times_VE11(dicom_list)
    VD13_prep_times = blood_tools.get_T2_prep_times_VD13A(dicom_list)

    if VE11_prep_times:
        prep_times = VE11_prep_times
    if VB17_prep_times:
        prep_times = VB17_prep_times
    if VD13_prep_times:
        prep_times = VD13_prep_times
        
     #for some reason, zero prep time image gets no prep time, assign it TE=2000    
    if len(prep_times) < len(dicom_list):
        prep_times=np.concatenate([np.array([2000]),prep_times])          

    signal = []
    log_signal = []
    std_signal = []
    slice_prep_times = []
    for prep_time, image, roi, slice_included in zip(prep_times, image_list, roi_list, included_slices):
        # handle slice exclusions
        if slice_included:
            slice_prep_times.append(prep_time)
            mean_signal_magnitude_over_ROI = blood_tools.calc_ROI_mean(
                roi, image)
            std_signal_magnitude_over_ROI = blood_tools.calc_ROI_std(
                roi, image)
            signal.append(mean_signal_magnitude_over_ROI)
            std_signal.append(std_signal_magnitude_over_ROI)
            log_signal.append(math.log(mean_signal_magnitude_over_ROI))
    if log_scale:
        return slice_prep_times, log_signal, std_signal
    else:
        return slice_prep_times, signal, std_signal


def get_T1_decay_signal(image_attributes, image_list, roi_list, included_slices, log_scale=False):
    """Gets T1 inversion times and signal"""
    inversion_times = [att['InversionTime'] for att in image_attributes]

    signal = []  # average signal across the ROI
    std = []
    log_signal = []  # log of the average signal
    slice_inversion_times = []
    for inv_time, image, roi, slice_included in zip(inversion_times, image_list, roi_list, included_slices):
        # handle slice exclusions
        if slice_included:
            mean_signal_magnitude_over_ROI = blood_tools.calc_ROI_mean(
                roi, image)
            std_signal = blood_tools.calc_ROI_std(roi, image)
            signal.append(mean_signal_magnitude_over_ROI)
            log_signal.append(math.log(mean_signal_magnitude_over_ROI))
            std.append(std_signal)
            slice_inversion_times.append(inv_time)

    if log_scale:
        return slice_inversion_times, np.array(log_signal), std
    else:
        return slice_inversion_times, np.array(signal), std


def bootstrap_decay_signal(image_attributes, image_list, roi_list, included_slices, iterations=1000):
    bootstrap_signals = []
    for ii in np.arange(iterations):
        sig = []
        for image, roi, slice_included in zip(image_list, roi_list, included_slices):
            if slice_included:
                pixels = roi.get_indices()
                npix = len(pixels[0])
                ind = np.random.randint(npix, size=npix)
                sig.append(image[pixels[0][ind], pixels[1][ind]].mean())
        bootstrap_signals.append(sig)
    return bootstrap_signals


def fit_bootstrap(x, sigs, model):
    # keep original best fit pars
    best_fit_pars = [p.value for p in model.pars]
    par_names = [p.name for p in model.pars]
    par_records = []
    for sig in sigs:
        try:
            model.fit(x, sig)
        except RuntimeError:
            continue
        par_records.append([p.value for p in model.pars])

        # reset back to best fit
        for qq, par in enumerate(model.pars):
            model[par.name].value = best_fit_pars[qq]  # reset to best fit pars

    par_uncertainties = {}
    for jj in np.arange(len(best_fit_pars)):
        par_hist = [record[jj] for record in par_records]
        par_uncertainties[par_names[jj]] = np.std(par_hist)

    return par_uncertainties


def main():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    win = MainWindow()
    win.show()
    app.exec_()
    return win


if __name__ == '__main__':
    win = main()
