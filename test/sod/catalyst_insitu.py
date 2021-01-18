
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

# Output frequency
output_freq = 500

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.8.1
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.8.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.8.1
      #
      # To ensure correct image size when batch processing, please search
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [2098, 1130]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [0.5, 0.1, 0.0]
      renderView1.StereoType = 'Crystal Eyes'
      renderView1.CameraPosition = [0.5, 0.1, 10000.0]
      renderView1.CameraFocalPoint = [0.5, 0.1, 0.0]
      renderView1.CameraFocalDisk = 1.0
      renderView1.CameraParallelScale = 0.5099019513592785

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='Density_%t.png', freq=output_freq, fittoscreen=0, magnification=1, width=2098, height=1130, cinema={}, compression=5)
      renderView1.ViewTime = datadescription.GetTime()

      SetActiveView(None)

      # ----------------------------------------------------------------
      # setup view layouts
      # ----------------------------------------------------------------

      # create new layout object 'Layout #1'
      layout1 = CreateLayout(name='Layout #1')
      layout1.AssignView(0, renderView1)

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Unstructured Grid Reader'
      # create a producer from a simulation input
      hydroLag = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from hydroLag
      hydroLagDisplay = Show(hydroLag, renderView1, 'UnstructuredGridRepresentation')

      # get color transfer function/color map for 'Density'
      densityLUT = GetColorTransferFunction('Density')
      densityLUT.RGBPoints = [0.125, 0.231373, 0.298039, 0.752941, 0.5625, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
      densityLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'Density'
      densityPWF = GetOpacityTransferFunction('Density')
      densityPWF.Points = [0.125, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
      densityPWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      hydroLagDisplay.Representation = 'Surface'
      hydroLagDisplay.ColorArrayName = ['CELLS', 'Density']
      hydroLagDisplay.LookupTable = densityLUT
      hydroLagDisplay.OSPRayScaleArray = 'NodeForce'
      hydroLagDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      hydroLagDisplay.SelectOrientationVectors = 'NodeForce'
      hydroLagDisplay.ScaleFactor = 0.1
      hydroLagDisplay.SelectScaleArray = 'NodeForce'
      hydroLagDisplay.GlyphType = 'Arrow'
      hydroLagDisplay.GlyphTableIndexArray = 'NodeForce'
      hydroLagDisplay.GaussianRadius = 0.005
      hydroLagDisplay.SetScaleArray = ['POINTS', 'NodeForce']
      hydroLagDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      hydroLagDisplay.OpacityArray = ['POINTS', 'NodeForce']
      hydroLagDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      hydroLagDisplay.DataAxesGrid = 'GridAxesRepresentation'
      hydroLagDisplay.PolarAxes = 'PolarAxesRepresentation'
      hydroLagDisplay.ScalarOpacityFunction = densityPWF
      hydroLagDisplay.ScalarOpacityUnitDistance = 0.12848724038000534

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      hydroLagDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      hydroLagDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for densityLUT in view renderView1
      densityLUTColorBar = GetScalarBar(densityLUT, renderView1)
      densityLUTColorBar.Title = 'Density'
      densityLUTColorBar.ComponentTitle = ''

      # set color bar visibility
      densityLUTColorBar.Visibility = 1

      # show color legend
      hydroLagDisplay.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(hydroLag)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [output_freq]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['NodeForce', 0], ['NodeMass', 0], ['NodeVelocity', 0], ['ArtificialViscosity', 1], ['Density', 1], ['InternalEnergy', 1], ['Mass', 1], ['Pressure', 1], ['SoundSpeed', 1], ['Volume', 1]]
    coprocessor.SetRequestedArrays('input', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
