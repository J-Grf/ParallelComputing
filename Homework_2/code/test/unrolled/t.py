# trace generated using paraview version 5.7.0-RC3
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
diskfinestvtu = XMLUnstructuredGridReader(FileName=['/hpcwork/pn744180/parallel_computing/Homework_2/code/test/unrolled/disk-finest.vtu'])
diskfinestvtu.CellArrayStatus = []
diskfinestvtu.PointArrayStatus = ['Temperature']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [895, 603]

# show data in view
diskfinestvtuDisplay = Show(diskfinestvtu, renderView1)

# get color transfer function/color map for 'Temperature'
temperatureLUT = GetColorTransferFunction('Temperature')
temperatureLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
temperatureLUT.InterpretValuesAsCategories = 0
temperatureLUT.AnnotationsInitialized = 0
temperatureLUT.ShowCategoricalColorsinDataRangeOnly = 0
temperatureLUT.RescaleOnVisibilityChange = 0
temperatureLUT.EnableOpacityMapping = 0
temperatureLUT.RGBPoints = [300.0, 0.231373, 0.298039, 0.752941, 302.21445272384386, 0.865003, 0.865003, 0.865003, 304.4289054476878, 0.705882, 0.0156863, 0.14902]
temperatureLUT.UseLogScale = 0
temperatureLUT.ColorSpace = 'Diverging'
temperatureLUT.UseBelowRangeColor = 0
temperatureLUT.BelowRangeColor = [0.0, 0.0, 0.0]
temperatureLUT.UseAboveRangeColor = 0
temperatureLUT.AboveRangeColor = [0.5, 0.5, 0.5]
temperatureLUT.NanColor = [1.0, 1.0, 0.0]
temperatureLUT.NanOpacity = 1.0
temperatureLUT.Discretize = 1
temperatureLUT.NumberOfTableValues = 256
temperatureLUT.ScalarRangeInitialized = 1.0
temperatureLUT.HSVWrap = 0
temperatureLUT.VectorComponent = 0
temperatureLUT.VectorMode = 'Magnitude'
temperatureLUT.AllowDuplicateScalars = 1
temperatureLUT.Annotations = []
temperatureLUT.ActiveAnnotatedValues = []
temperatureLUT.IndexedColors = []
temperatureLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Temperature'
temperaturePWF = GetOpacityTransferFunction('Temperature')
temperaturePWF.Points = [300.0, 0.0, 0.5, 0.0, 304.4289054476878, 1.0, 0.5, 0.0]
temperaturePWF.AllowDuplicateScalars = 1
temperaturePWF.UseLogScale = 0
temperaturePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
diskfinestvtuDisplay.Representation = 'Surface'
diskfinestvtuDisplay.AmbientColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.ColorArrayName = ['POINTS', 'Temperature']
diskfinestvtuDisplay.DiffuseColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.LookupTable = temperatureLUT
diskfinestvtuDisplay.MapScalars = 1
diskfinestvtuDisplay.MultiComponentsMapping = 0
diskfinestvtuDisplay.InterpolateScalarsBeforeMapping = 1
diskfinestvtuDisplay.Opacity = 1.0
diskfinestvtuDisplay.PointSize = 2.0
diskfinestvtuDisplay.LineWidth = 1.0
diskfinestvtuDisplay.RenderLinesAsTubes = 0
diskfinestvtuDisplay.RenderPointsAsSpheres = 0
diskfinestvtuDisplay.Interpolation = 'Gouraud'
diskfinestvtuDisplay.Specular = 0.0
diskfinestvtuDisplay.SpecularColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.SpecularPower = 100.0
diskfinestvtuDisplay.Luminosity = 0.0
diskfinestvtuDisplay.Ambient = 0.0
diskfinestvtuDisplay.Diffuse = 1.0
diskfinestvtuDisplay.EdgeColor = [0.0, 0.0, 0.5]
diskfinestvtuDisplay.BackfaceRepresentation = 'Follow Frontface'
diskfinestvtuDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.BackfaceOpacity = 1.0
diskfinestvtuDisplay.Position = [0.0, 0.0, 0.0]
diskfinestvtuDisplay.Scale = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.Orientation = [0.0, 0.0, 0.0]
diskfinestvtuDisplay.Origin = [0.0, 0.0, 0.0]
diskfinestvtuDisplay.Pickable = 1
diskfinestvtuDisplay.Texture = None
diskfinestvtuDisplay.Triangulate = 0
diskfinestvtuDisplay.UseShaderReplacements = 0
diskfinestvtuDisplay.ShaderReplacements = ''
diskfinestvtuDisplay.NonlinearSubdivisionLevel = 1
diskfinestvtuDisplay.UseDataPartitions = 0
diskfinestvtuDisplay.OSPRayUseScaleArray = 0
diskfinestvtuDisplay.OSPRayScaleArray = 'Temperature'
diskfinestvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
diskfinestvtuDisplay.OSPRayMaterial = 'None'
diskfinestvtuDisplay.Orient = 0
diskfinestvtuDisplay.OrientationMode = 'Direction'
diskfinestvtuDisplay.SelectOrientationVectors = 'None'
diskfinestvtuDisplay.Scaling = 0
diskfinestvtuDisplay.ScaleMode = 'No Data Scaling Off'
diskfinestvtuDisplay.ScaleFactor = 0.020000000000000004
diskfinestvtuDisplay.SelectScaleArray = 'Temperature'
diskfinestvtuDisplay.GlyphType = 'Arrow'
diskfinestvtuDisplay.UseGlyphTable = 0
diskfinestvtuDisplay.GlyphTableIndexArray = 'Temperature'
diskfinestvtuDisplay.UseCompositeGlyphTable = 0
diskfinestvtuDisplay.UseGlyphCullingAndLOD = 0
diskfinestvtuDisplay.LODValues = []
diskfinestvtuDisplay.ColorByLODIndex = 0
diskfinestvtuDisplay.GaussianRadius = 0.001
diskfinestvtuDisplay.ShaderPreset = 'Sphere'
diskfinestvtuDisplay.CustomTriangleScale = 3
diskfinestvtuDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
diskfinestvtuDisplay.Emissive = 0
diskfinestvtuDisplay.ScaleByArray = 0
diskfinestvtuDisplay.SetScaleArray = ['POINTS', 'Temperature']
diskfinestvtuDisplay.ScaleArrayComponent = ''
diskfinestvtuDisplay.UseScaleFunction = 1
diskfinestvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
diskfinestvtuDisplay.OpacityByArray = 0
diskfinestvtuDisplay.OpacityArray = ['POINTS', 'Temperature']
diskfinestvtuDisplay.OpacityArrayComponent = ''
diskfinestvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
diskfinestvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
diskfinestvtuDisplay.SelectionCellLabelBold = 0
diskfinestvtuDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
diskfinestvtuDisplay.SelectionCellLabelFontFamily = 'Arial'
diskfinestvtuDisplay.SelectionCellLabelFontFile = ''
diskfinestvtuDisplay.SelectionCellLabelFontSize = 18
diskfinestvtuDisplay.SelectionCellLabelItalic = 0
diskfinestvtuDisplay.SelectionCellLabelJustification = 'Left'
diskfinestvtuDisplay.SelectionCellLabelOpacity = 1.0
diskfinestvtuDisplay.SelectionCellLabelShadow = 0
diskfinestvtuDisplay.SelectionPointLabelBold = 0
diskfinestvtuDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
diskfinestvtuDisplay.SelectionPointLabelFontFamily = 'Arial'
diskfinestvtuDisplay.SelectionPointLabelFontFile = ''
diskfinestvtuDisplay.SelectionPointLabelFontSize = 18
diskfinestvtuDisplay.SelectionPointLabelItalic = 0
diskfinestvtuDisplay.SelectionPointLabelJustification = 'Left'
diskfinestvtuDisplay.SelectionPointLabelOpacity = 1.0
diskfinestvtuDisplay.SelectionPointLabelShadow = 0
diskfinestvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
diskfinestvtuDisplay.ScalarOpacityFunction = temperaturePWF
diskfinestvtuDisplay.ScalarOpacityUnitDistance = 0.006441580108262787
diskfinestvtuDisplay.ExtractedBlockIndex = 0
diskfinestvtuDisplay.SelectMapper = 'Projected tetra'
diskfinestvtuDisplay.SamplingDimensions = [128, 128, 128]
diskfinestvtuDisplay.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
diskfinestvtuDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
diskfinestvtuDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
diskfinestvtuDisplay.GlyphType.TipResolution = 6
diskfinestvtuDisplay.GlyphType.TipRadius = 0.1
diskfinestvtuDisplay.GlyphType.TipLength = 0.35
diskfinestvtuDisplay.GlyphType.ShaftResolution = 6
diskfinestvtuDisplay.GlyphType.ShaftRadius = 0.03
diskfinestvtuDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
diskfinestvtuDisplay.ScaleTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 304.4289054476878, 1.0, 0.5, 0.0]
diskfinestvtuDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
diskfinestvtuDisplay.OpacityTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 304.4289054476878, 1.0, 0.5, 0.0]
diskfinestvtuDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
diskfinestvtuDisplay.DataAxesGrid.XTitle = 'X Axis'
diskfinestvtuDisplay.DataAxesGrid.YTitle = 'Y Axis'
diskfinestvtuDisplay.DataAxesGrid.ZTitle = 'Z Axis'
diskfinestvtuDisplay.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
diskfinestvtuDisplay.DataAxesGrid.XTitleFontFile = ''
diskfinestvtuDisplay.DataAxesGrid.XTitleBold = 0
diskfinestvtuDisplay.DataAxesGrid.XTitleItalic = 0
diskfinestvtuDisplay.DataAxesGrid.XTitleFontSize = 12
diskfinestvtuDisplay.DataAxesGrid.XTitleShadow = 0
diskfinestvtuDisplay.DataAxesGrid.XTitleOpacity = 1.0
diskfinestvtuDisplay.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
diskfinestvtuDisplay.DataAxesGrid.YTitleFontFile = ''
diskfinestvtuDisplay.DataAxesGrid.YTitleBold = 0
diskfinestvtuDisplay.DataAxesGrid.YTitleItalic = 0
diskfinestvtuDisplay.DataAxesGrid.YTitleFontSize = 12
diskfinestvtuDisplay.DataAxesGrid.YTitleShadow = 0
diskfinestvtuDisplay.DataAxesGrid.YTitleOpacity = 1.0
diskfinestvtuDisplay.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
diskfinestvtuDisplay.DataAxesGrid.ZTitleFontFile = ''
diskfinestvtuDisplay.DataAxesGrid.ZTitleBold = 0
diskfinestvtuDisplay.DataAxesGrid.ZTitleItalic = 0
diskfinestvtuDisplay.DataAxesGrid.ZTitleFontSize = 12
diskfinestvtuDisplay.DataAxesGrid.ZTitleShadow = 0
diskfinestvtuDisplay.DataAxesGrid.ZTitleOpacity = 1.0
diskfinestvtuDisplay.DataAxesGrid.FacesToRender = 63
diskfinestvtuDisplay.DataAxesGrid.CullBackface = 0
diskfinestvtuDisplay.DataAxesGrid.CullFrontface = 1
diskfinestvtuDisplay.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.DataAxesGrid.ShowGrid = 0
diskfinestvtuDisplay.DataAxesGrid.ShowEdges = 1
diskfinestvtuDisplay.DataAxesGrid.ShowTicks = 1
diskfinestvtuDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
diskfinestvtuDisplay.DataAxesGrid.AxesToLabel = 63
diskfinestvtuDisplay.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
diskfinestvtuDisplay.DataAxesGrid.XLabelFontFile = ''
diskfinestvtuDisplay.DataAxesGrid.XLabelBold = 0
diskfinestvtuDisplay.DataAxesGrid.XLabelItalic = 0
diskfinestvtuDisplay.DataAxesGrid.XLabelFontSize = 12
diskfinestvtuDisplay.DataAxesGrid.XLabelShadow = 0
diskfinestvtuDisplay.DataAxesGrid.XLabelOpacity = 1.0
diskfinestvtuDisplay.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
diskfinestvtuDisplay.DataAxesGrid.YLabelFontFile = ''
diskfinestvtuDisplay.DataAxesGrid.YLabelBold = 0
diskfinestvtuDisplay.DataAxesGrid.YLabelItalic = 0
diskfinestvtuDisplay.DataAxesGrid.YLabelFontSize = 12
diskfinestvtuDisplay.DataAxesGrid.YLabelShadow = 0
diskfinestvtuDisplay.DataAxesGrid.YLabelOpacity = 1.0
diskfinestvtuDisplay.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
diskfinestvtuDisplay.DataAxesGrid.ZLabelFontFile = ''
diskfinestvtuDisplay.DataAxesGrid.ZLabelBold = 0
diskfinestvtuDisplay.DataAxesGrid.ZLabelItalic = 0
diskfinestvtuDisplay.DataAxesGrid.ZLabelFontSize = 12
diskfinestvtuDisplay.DataAxesGrid.ZLabelShadow = 0
diskfinestvtuDisplay.DataAxesGrid.ZLabelOpacity = 1.0
diskfinestvtuDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
diskfinestvtuDisplay.DataAxesGrid.XAxisPrecision = 2
diskfinestvtuDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
diskfinestvtuDisplay.DataAxesGrid.XAxisLabels = []
diskfinestvtuDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
diskfinestvtuDisplay.DataAxesGrid.YAxisPrecision = 2
diskfinestvtuDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
diskfinestvtuDisplay.DataAxesGrid.YAxisLabels = []
diskfinestvtuDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
diskfinestvtuDisplay.DataAxesGrid.ZAxisPrecision = 2
diskfinestvtuDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
diskfinestvtuDisplay.DataAxesGrid.ZAxisLabels = []
diskfinestvtuDisplay.DataAxesGrid.UseCustomBounds = 0
diskfinestvtuDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
diskfinestvtuDisplay.PolarAxes.Visibility = 0
diskfinestvtuDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
diskfinestvtuDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
diskfinestvtuDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
diskfinestvtuDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
diskfinestvtuDisplay.PolarAxes.EnableCustomRange = 0
diskfinestvtuDisplay.PolarAxes.CustomRange = [0.0, 1.0]
diskfinestvtuDisplay.PolarAxes.PolarAxisVisibility = 1
diskfinestvtuDisplay.PolarAxes.RadialAxesVisibility = 1
diskfinestvtuDisplay.PolarAxes.DrawRadialGridlines = 1
diskfinestvtuDisplay.PolarAxes.PolarArcsVisibility = 1
diskfinestvtuDisplay.PolarAxes.DrawPolarArcsGridlines = 1
diskfinestvtuDisplay.PolarAxes.NumberOfRadialAxes = 0
diskfinestvtuDisplay.PolarAxes.AutoSubdividePolarAxis = 1
diskfinestvtuDisplay.PolarAxes.NumberOfPolarAxis = 0
diskfinestvtuDisplay.PolarAxes.MinimumRadius = 0.0
diskfinestvtuDisplay.PolarAxes.MinimumAngle = 0.0
diskfinestvtuDisplay.PolarAxes.MaximumAngle = 90.0
diskfinestvtuDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
diskfinestvtuDisplay.PolarAxes.Ratio = 1.0
diskfinestvtuDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleVisibility = 1
diskfinestvtuDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
diskfinestvtuDisplay.PolarAxes.PolarLabelVisibility = 1
diskfinestvtuDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
diskfinestvtuDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
diskfinestvtuDisplay.PolarAxes.RadialLabelVisibility = 1
diskfinestvtuDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
diskfinestvtuDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
diskfinestvtuDisplay.PolarAxes.RadialUnitsVisibility = 1
diskfinestvtuDisplay.PolarAxes.ScreenSize = 10.0
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleFontFile = ''
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleBold = 0
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleItalic = 0
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleShadow = 0
diskfinestvtuDisplay.PolarAxes.PolarAxisTitleFontSize = 12
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelFontFile = ''
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelBold = 0
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelItalic = 0
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelShadow = 0
diskfinestvtuDisplay.PolarAxes.PolarAxisLabelFontSize = 12
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextBold = 0
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextItalic = 0
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextShadow = 0
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
diskfinestvtuDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
diskfinestvtuDisplay.PolarAxes.EnableDistanceLOD = 1
diskfinestvtuDisplay.PolarAxes.DistanceLODThreshold = 0.7
diskfinestvtuDisplay.PolarAxes.EnableViewAngleLOD = 1
diskfinestvtuDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
diskfinestvtuDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
diskfinestvtuDisplay.PolarAxes.PolarTicksVisibility = 1
diskfinestvtuDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
diskfinestvtuDisplay.PolarAxes.TickLocation = 'Both'
diskfinestvtuDisplay.PolarAxes.AxisTickVisibility = 1
diskfinestvtuDisplay.PolarAxes.AxisMinorTickVisibility = 0
diskfinestvtuDisplay.PolarAxes.ArcTickVisibility = 1
diskfinestvtuDisplay.PolarAxes.ArcMinorTickVisibility = 0
diskfinestvtuDisplay.PolarAxes.DeltaAngleMajor = 10.0
diskfinestvtuDisplay.PolarAxes.DeltaAngleMinor = 5.0
diskfinestvtuDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
diskfinestvtuDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
diskfinestvtuDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
diskfinestvtuDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
diskfinestvtuDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
diskfinestvtuDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
diskfinestvtuDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
diskfinestvtuDisplay.PolarAxes.ArcMajorTickSize = 0.0
diskfinestvtuDisplay.PolarAxes.ArcTickRatioSize = 0.3
diskfinestvtuDisplay.PolarAxes.ArcMajorTickThickness = 1.0
diskfinestvtuDisplay.PolarAxes.ArcTickRatioThickness = 0.5
diskfinestvtuDisplay.PolarAxes.Use2DMode = 0
diskfinestvtuDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.0, 0.0, 10000.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
diskfinestvtuDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=diskfinestvtu,
    Source='High Resolution Line Source')
plotOverLine1.PassPartialArrays = 1
plotOverLine1.ComputeTolerance = 1
plotOverLine1.Tolerance = 2.220446049250313e-16

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine1.Source.Point1 = [-0.1, -0.1, 0.0]
plotOverLine1.Source.Point2 = [0.1, 0.1, 0.0]
plotOverLine1.Source.Resolution = 1000

# show data in view
plotOverLine1Display = Show(plotOverLine1, renderView1)

# trace defaults for the display properties.
plotOverLine1Display.Representation = 'Surface'
plotOverLine1Display.AmbientColor = [1.0, 1.0, 1.0]
plotOverLine1Display.ColorArrayName = ['POINTS', 'Temperature']
plotOverLine1Display.DiffuseColor = [1.0, 1.0, 1.0]
plotOverLine1Display.LookupTable = temperatureLUT
plotOverLine1Display.MapScalars = 1
plotOverLine1Display.MultiComponentsMapping = 0
plotOverLine1Display.InterpolateScalarsBeforeMapping = 1
plotOverLine1Display.Opacity = 1.0
plotOverLine1Display.PointSize = 2.0
plotOverLine1Display.LineWidth = 1.0
plotOverLine1Display.RenderLinesAsTubes = 0
plotOverLine1Display.RenderPointsAsSpheres = 0
plotOverLine1Display.Interpolation = 'Gouraud'
plotOverLine1Display.Specular = 0.0
plotOverLine1Display.SpecularColor = [1.0, 1.0, 1.0]
plotOverLine1Display.SpecularPower = 100.0
plotOverLine1Display.Luminosity = 0.0
plotOverLine1Display.Ambient = 0.0
plotOverLine1Display.Diffuse = 1.0
plotOverLine1Display.EdgeColor = [0.0, 0.0, 0.5]
plotOverLine1Display.BackfaceRepresentation = 'Follow Frontface'
plotOverLine1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
plotOverLine1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
plotOverLine1Display.BackfaceOpacity = 1.0
plotOverLine1Display.Position = [0.0, 0.0, 0.0]
plotOverLine1Display.Scale = [1.0, 1.0, 1.0]
plotOverLine1Display.Orientation = [0.0, 0.0, 0.0]
plotOverLine1Display.Origin = [0.0, 0.0, 0.0]
plotOverLine1Display.Pickable = 1
plotOverLine1Display.Texture = None
plotOverLine1Display.Triangulate = 0
plotOverLine1Display.UseShaderReplacements = 0
plotOverLine1Display.ShaderReplacements = ''
plotOverLine1Display.NonlinearSubdivisionLevel = 1
plotOverLine1Display.UseDataPartitions = 0
plotOverLine1Display.OSPRayUseScaleArray = 0
plotOverLine1Display.OSPRayScaleArray = 'Temperature'
plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plotOverLine1Display.OSPRayMaterial = 'None'
plotOverLine1Display.Orient = 0
plotOverLine1Display.OrientationMode = 'Direction'
plotOverLine1Display.SelectOrientationVectors = 'None'
plotOverLine1Display.Scaling = 0
plotOverLine1Display.ScaleMode = 'No Data Scaling Off'
plotOverLine1Display.ScaleFactor = 0.020000000298023225
plotOverLine1Display.SelectScaleArray = 'Temperature'
plotOverLine1Display.GlyphType = 'Arrow'
plotOverLine1Display.UseGlyphTable = 0
plotOverLine1Display.GlyphTableIndexArray = 'Temperature'
plotOverLine1Display.UseCompositeGlyphTable = 0
plotOverLine1Display.UseGlyphCullingAndLOD = 0
plotOverLine1Display.LODValues = []
plotOverLine1Display.ColorByLODIndex = 0
plotOverLine1Display.GaussianRadius = 0.0010000000149011613
plotOverLine1Display.ShaderPreset = 'Sphere'
plotOverLine1Display.CustomTriangleScale = 3
plotOverLine1Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
plotOverLine1Display.Emissive = 0
plotOverLine1Display.ScaleByArray = 0
plotOverLine1Display.SetScaleArray = ['POINTS', 'Temperature']
plotOverLine1Display.ScaleArrayComponent = ''
plotOverLine1Display.UseScaleFunction = 1
plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.OpacityByArray = 0
plotOverLine1Display.OpacityArray = ['POINTS', 'Temperature']
plotOverLine1Display.OpacityArrayComponent = ''
plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
plotOverLine1Display.SelectionCellLabelBold = 0
plotOverLine1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
plotOverLine1Display.SelectionCellLabelFontFamily = 'Arial'
plotOverLine1Display.SelectionCellLabelFontFile = ''
plotOverLine1Display.SelectionCellLabelFontSize = 18
plotOverLine1Display.SelectionCellLabelItalic = 0
plotOverLine1Display.SelectionCellLabelJustification = 'Left'
plotOverLine1Display.SelectionCellLabelOpacity = 1.0
plotOverLine1Display.SelectionCellLabelShadow = 0
plotOverLine1Display.SelectionPointLabelBold = 0
plotOverLine1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
plotOverLine1Display.SelectionPointLabelFontFamily = 'Arial'
plotOverLine1Display.SelectionPointLabelFontFile = ''
plotOverLine1Display.SelectionPointLabelFontSize = 18
plotOverLine1Display.SelectionPointLabelItalic = 0
plotOverLine1Display.SelectionPointLabelJustification = 'Left'
plotOverLine1Display.SelectionPointLabelOpacity = 1.0
plotOverLine1Display.SelectionPointLabelShadow = 0
plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
plotOverLine1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
plotOverLine1Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
plotOverLine1Display.GlyphType.TipResolution = 6
plotOverLine1Display.GlyphType.TipRadius = 0.1
plotOverLine1Display.GlyphType.TipLength = 0.35
plotOverLine1Display.GlyphType.ShaftResolution = 6
plotOverLine1Display.GlyphType.ShaftRadius = 0.03
plotOverLine1Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plotOverLine1Display.ScaleTransferFunction.Points = [300.0024434411071, 0.0, 0.5, 0.0, 304.4289054476878, 1.0, 0.5, 0.0]
plotOverLine1Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plotOverLine1Display.OpacityTransferFunction.Points = [300.0024434411071, 0.0, 0.5, 0.0, 304.4289054476878, 1.0, 0.5, 0.0]
plotOverLine1Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
plotOverLine1Display.DataAxesGrid.XTitle = 'X Axis'
plotOverLine1Display.DataAxesGrid.YTitle = 'Y Axis'
plotOverLine1Display.DataAxesGrid.ZTitle = 'Z Axis'
plotOverLine1Display.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
plotOverLine1Display.DataAxesGrid.XTitleFontFamily = 'Arial'
plotOverLine1Display.DataAxesGrid.XTitleFontFile = ''
plotOverLine1Display.DataAxesGrid.XTitleBold = 0
plotOverLine1Display.DataAxesGrid.XTitleItalic = 0
plotOverLine1Display.DataAxesGrid.XTitleFontSize = 12
plotOverLine1Display.DataAxesGrid.XTitleShadow = 0
plotOverLine1Display.DataAxesGrid.XTitleOpacity = 1.0
plotOverLine1Display.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
plotOverLine1Display.DataAxesGrid.YTitleFontFamily = 'Arial'
plotOverLine1Display.DataAxesGrid.YTitleFontFile = ''
plotOverLine1Display.DataAxesGrid.YTitleBold = 0
plotOverLine1Display.DataAxesGrid.YTitleItalic = 0
plotOverLine1Display.DataAxesGrid.YTitleFontSize = 12
plotOverLine1Display.DataAxesGrid.YTitleShadow = 0
plotOverLine1Display.DataAxesGrid.YTitleOpacity = 1.0
plotOverLine1Display.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
plotOverLine1Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
plotOverLine1Display.DataAxesGrid.ZTitleFontFile = ''
plotOverLine1Display.DataAxesGrid.ZTitleBold = 0
plotOverLine1Display.DataAxesGrid.ZTitleItalic = 0
plotOverLine1Display.DataAxesGrid.ZTitleFontSize = 12
plotOverLine1Display.DataAxesGrid.ZTitleShadow = 0
plotOverLine1Display.DataAxesGrid.ZTitleOpacity = 1.0
plotOverLine1Display.DataAxesGrid.FacesToRender = 63
plotOverLine1Display.DataAxesGrid.CullBackface = 0
plotOverLine1Display.DataAxesGrid.CullFrontface = 1
plotOverLine1Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
plotOverLine1Display.DataAxesGrid.ShowGrid = 0
plotOverLine1Display.DataAxesGrid.ShowEdges = 1
plotOverLine1Display.DataAxesGrid.ShowTicks = 1
plotOverLine1Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
plotOverLine1Display.DataAxesGrid.AxesToLabel = 63
plotOverLine1Display.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
plotOverLine1Display.DataAxesGrid.XLabelFontFamily = 'Arial'
plotOverLine1Display.DataAxesGrid.XLabelFontFile = ''
plotOverLine1Display.DataAxesGrid.XLabelBold = 0
plotOverLine1Display.DataAxesGrid.XLabelItalic = 0
plotOverLine1Display.DataAxesGrid.XLabelFontSize = 12
plotOverLine1Display.DataAxesGrid.XLabelShadow = 0
plotOverLine1Display.DataAxesGrid.XLabelOpacity = 1.0
plotOverLine1Display.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
plotOverLine1Display.DataAxesGrid.YLabelFontFamily = 'Arial'
plotOverLine1Display.DataAxesGrid.YLabelFontFile = ''
plotOverLine1Display.DataAxesGrid.YLabelBold = 0
plotOverLine1Display.DataAxesGrid.YLabelItalic = 0
plotOverLine1Display.DataAxesGrid.YLabelFontSize = 12
plotOverLine1Display.DataAxesGrid.YLabelShadow = 0
plotOverLine1Display.DataAxesGrid.YLabelOpacity = 1.0
plotOverLine1Display.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
plotOverLine1Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
plotOverLine1Display.DataAxesGrid.ZLabelFontFile = ''
plotOverLine1Display.DataAxesGrid.ZLabelBold = 0
plotOverLine1Display.DataAxesGrid.ZLabelItalic = 0
plotOverLine1Display.DataAxesGrid.ZLabelFontSize = 12
plotOverLine1Display.DataAxesGrid.ZLabelShadow = 0
plotOverLine1Display.DataAxesGrid.ZLabelOpacity = 1.0
plotOverLine1Display.DataAxesGrid.XAxisNotation = 'Mixed'
plotOverLine1Display.DataAxesGrid.XAxisPrecision = 2
plotOverLine1Display.DataAxesGrid.XAxisUseCustomLabels = 0
plotOverLine1Display.DataAxesGrid.XAxisLabels = []
plotOverLine1Display.DataAxesGrid.YAxisNotation = 'Mixed'
plotOverLine1Display.DataAxesGrid.YAxisPrecision = 2
plotOverLine1Display.DataAxesGrid.YAxisUseCustomLabels = 0
plotOverLine1Display.DataAxesGrid.YAxisLabels = []
plotOverLine1Display.DataAxesGrid.ZAxisNotation = 'Mixed'
plotOverLine1Display.DataAxesGrid.ZAxisPrecision = 2
plotOverLine1Display.DataAxesGrid.ZAxisUseCustomLabels = 0
plotOverLine1Display.DataAxesGrid.ZAxisLabels = []
plotOverLine1Display.DataAxesGrid.UseCustomBounds = 0
plotOverLine1Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
plotOverLine1Display.PolarAxes.Visibility = 0
plotOverLine1Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
plotOverLine1Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
plotOverLine1Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
plotOverLine1Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
plotOverLine1Display.PolarAxes.EnableCustomRange = 0
plotOverLine1Display.PolarAxes.CustomRange = [0.0, 1.0]
plotOverLine1Display.PolarAxes.PolarAxisVisibility = 1
plotOverLine1Display.PolarAxes.RadialAxesVisibility = 1
plotOverLine1Display.PolarAxes.DrawRadialGridlines = 1
plotOverLine1Display.PolarAxes.PolarArcsVisibility = 1
plotOverLine1Display.PolarAxes.DrawPolarArcsGridlines = 1
plotOverLine1Display.PolarAxes.NumberOfRadialAxes = 0
plotOverLine1Display.PolarAxes.AutoSubdividePolarAxis = 1
plotOverLine1Display.PolarAxes.NumberOfPolarAxis = 0
plotOverLine1Display.PolarAxes.MinimumRadius = 0.0
plotOverLine1Display.PolarAxes.MinimumAngle = 0.0
plotOverLine1Display.PolarAxes.MaximumAngle = 90.0
plotOverLine1Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
plotOverLine1Display.PolarAxes.Ratio = 1.0
plotOverLine1Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.PolarAxisTitleVisibility = 1
plotOverLine1Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
plotOverLine1Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
plotOverLine1Display.PolarAxes.PolarLabelVisibility = 1
plotOverLine1Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
plotOverLine1Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
plotOverLine1Display.PolarAxes.RadialLabelVisibility = 1
plotOverLine1Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
plotOverLine1Display.PolarAxes.RadialLabelLocation = 'Bottom'
plotOverLine1Display.PolarAxes.RadialUnitsVisibility = 1
plotOverLine1Display.PolarAxes.ScreenSize = 10.0
plotOverLine1Display.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.PolarAxisTitleOpacity = 1.0
plotOverLine1Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
plotOverLine1Display.PolarAxes.PolarAxisTitleFontFile = ''
plotOverLine1Display.PolarAxes.PolarAxisTitleBold = 0
plotOverLine1Display.PolarAxes.PolarAxisTitleItalic = 0
plotOverLine1Display.PolarAxes.PolarAxisTitleShadow = 0
plotOverLine1Display.PolarAxes.PolarAxisTitleFontSize = 12
plotOverLine1Display.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.PolarAxisLabelOpacity = 1.0
plotOverLine1Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
plotOverLine1Display.PolarAxes.PolarAxisLabelFontFile = ''
plotOverLine1Display.PolarAxes.PolarAxisLabelBold = 0
plotOverLine1Display.PolarAxes.PolarAxisLabelItalic = 0
plotOverLine1Display.PolarAxes.PolarAxisLabelShadow = 0
plotOverLine1Display.PolarAxes.PolarAxisLabelFontSize = 12
plotOverLine1Display.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
plotOverLine1Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
plotOverLine1Display.PolarAxes.LastRadialAxisTextFontFile = ''
plotOverLine1Display.PolarAxes.LastRadialAxisTextBold = 0
plotOverLine1Display.PolarAxes.LastRadialAxisTextItalic = 0
plotOverLine1Display.PolarAxes.LastRadialAxisTextShadow = 0
plotOverLine1Display.PolarAxes.LastRadialAxisTextFontSize = 12
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextBold = 0
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
plotOverLine1Display.PolarAxes.EnableDistanceLOD = 1
plotOverLine1Display.PolarAxes.DistanceLODThreshold = 0.7
plotOverLine1Display.PolarAxes.EnableViewAngleLOD = 1
plotOverLine1Display.PolarAxes.ViewAngleLODThreshold = 0.7
plotOverLine1Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
plotOverLine1Display.PolarAxes.PolarTicksVisibility = 1
plotOverLine1Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
plotOverLine1Display.PolarAxes.TickLocation = 'Both'
plotOverLine1Display.PolarAxes.AxisTickVisibility = 1
plotOverLine1Display.PolarAxes.AxisMinorTickVisibility = 0
plotOverLine1Display.PolarAxes.ArcTickVisibility = 1
plotOverLine1Display.PolarAxes.ArcMinorTickVisibility = 0
plotOverLine1Display.PolarAxes.DeltaAngleMajor = 10.0
plotOverLine1Display.PolarAxes.DeltaAngleMinor = 5.0
plotOverLine1Display.PolarAxes.PolarAxisMajorTickSize = 0.0
plotOverLine1Display.PolarAxes.PolarAxisTickRatioSize = 0.3
plotOverLine1Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
plotOverLine1Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
plotOverLine1Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
plotOverLine1Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
plotOverLine1Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
plotOverLine1Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
plotOverLine1Display.PolarAxes.ArcMajorTickSize = 0.0
plotOverLine1Display.PolarAxes.ArcTickRatioSize = 0.3
plotOverLine1Display.PolarAxes.ArcMajorTickThickness = 1.0
plotOverLine1Display.PolarAxes.ArcTickRatioThickness = 0.5
plotOverLine1Display.PolarAxes.Use2DMode = 0
plotOverLine1Display.PolarAxes.UseLogAxis = 0

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.UseCache = 0
lineChartView1.ViewSize = [400, 400]
lineChartView1.ChartTitle = ''
lineChartView1.ChartTitleAlignment = 'Center'
lineChartView1.ChartTitleFontFamily = 'Arial'
lineChartView1.ChartTitleFontFile = ''
lineChartView1.ChartTitleFontSize = 18
lineChartView1.ChartTitleBold = 0
lineChartView1.ChartTitleItalic = 0
lineChartView1.ChartTitleColor = [0.0, 0.0, 0.0]
lineChartView1.ShowLegend = 1
lineChartView1.LegendLocation = 'TopRight'
lineChartView1.SortByXAxis = 0
lineChartView1.LegendPosition = [0, 0]
lineChartView1.LegendFontFamily = 'Arial'
lineChartView1.LegendFontFile = ''
lineChartView1.LegendFontSize = 12
lineChartView1.LegendBold = 0
lineChartView1.LegendItalic = 0
lineChartView1.TooltipNotation = 'Mixed'
lineChartView1.TooltipPrecision = 6
lineChartView1.HideTimeMarker = 0
lineChartView1.LeftAxisTitle = ''
lineChartView1.ShowLeftAxisGrid = 1
lineChartView1.LeftAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.LeftAxisColor = [0.0, 0.0, 0.0]
lineChartView1.LeftAxisTitleFontFamily = 'Arial'
lineChartView1.LeftAxisTitleFontFile = ''
lineChartView1.LeftAxisTitleFontSize = 18
lineChartView1.LeftAxisTitleBold = 1
lineChartView1.LeftAxisTitleItalic = 0
lineChartView1.LeftAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.LeftAxisLogScale = 0
lineChartView1.LeftAxisUseCustomRange = 0
lineChartView1.LeftAxisRangeMinimum = 0.0
lineChartView1.LeftAxisRangeMaximum = 1.0
lineChartView1.ShowLeftAxisLabels = 1
lineChartView1.LeftAxisLabelNotation = 'Mixed'
lineChartView1.LeftAxisLabelPrecision = 2
lineChartView1.LeftAxisUseCustomLabels = 0
lineChartView1.LeftAxisLabels = []
lineChartView1.LeftAxisLabelFontFamily = 'Arial'
lineChartView1.LeftAxisLabelFontFile = ''
lineChartView1.LeftAxisLabelFontSize = 12
lineChartView1.LeftAxisLabelBold = 0
lineChartView1.LeftAxisLabelItalic = 0
lineChartView1.LeftAxisLabelColor = [0.0, 0.0, 0.0]
lineChartView1.BottomAxisTitle = ''
lineChartView1.ShowBottomAxisGrid = 1
lineChartView1.BottomAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.BottomAxisColor = [0.0, 0.0, 0.0]
lineChartView1.BottomAxisTitleFontFamily = 'Arial'
lineChartView1.BottomAxisTitleFontFile = ''
lineChartView1.BottomAxisTitleFontSize = 18
lineChartView1.BottomAxisTitleBold = 1
lineChartView1.BottomAxisTitleItalic = 0
lineChartView1.BottomAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.BottomAxisLogScale = 0
lineChartView1.BottomAxisUseCustomRange = 0
lineChartView1.BottomAxisRangeMinimum = 0.0
lineChartView1.BottomAxisRangeMaximum = 1.0
lineChartView1.ShowBottomAxisLabels = 1
lineChartView1.BottomAxisLabelNotation = 'Mixed'
lineChartView1.BottomAxisLabelPrecision = 2
lineChartView1.BottomAxisUseCustomLabels = 0
lineChartView1.BottomAxisLabels = []
lineChartView1.BottomAxisLabelFontFamily = 'Arial'
lineChartView1.BottomAxisLabelFontFile = ''
lineChartView1.BottomAxisLabelFontSize = 12
lineChartView1.BottomAxisLabelBold = 0
lineChartView1.BottomAxisLabelItalic = 0
lineChartView1.BottomAxisLabelColor = [0.0, 0.0, 0.0]
lineChartView1.RightAxisTitle = ''
lineChartView1.ShowRightAxisGrid = 1
lineChartView1.RightAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.RightAxisColor = [0.0, 0.0, 0.0]
lineChartView1.RightAxisTitleFontFamily = 'Arial'
lineChartView1.RightAxisTitleFontFile = ''
lineChartView1.RightAxisTitleFontSize = 18
lineChartView1.RightAxisTitleBold = 1
lineChartView1.RightAxisTitleItalic = 0
lineChartView1.RightAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.RightAxisLogScale = 0
lineChartView1.RightAxisUseCustomRange = 0
lineChartView1.RightAxisRangeMinimum = 0.0
lineChartView1.RightAxisRangeMaximum = 1.0
lineChartView1.ShowRightAxisLabels = 1
lineChartView1.RightAxisLabelNotation = 'Mixed'
lineChartView1.RightAxisLabelPrecision = 2
lineChartView1.RightAxisUseCustomLabels = 0
lineChartView1.RightAxisLabels = []
lineChartView1.RightAxisLabelFontFamily = 'Arial'
lineChartView1.RightAxisLabelFontFile = ''
lineChartView1.RightAxisLabelFontSize = 12
lineChartView1.RightAxisLabelBold = 0
lineChartView1.RightAxisLabelItalic = 0
lineChartView1.RightAxisLabelColor = [0.0, 0.0, 0.0]
lineChartView1.TopAxisTitle = ''
lineChartView1.ShowTopAxisGrid = 1
lineChartView1.TopAxisGridColor = [0.95, 0.95, 0.95]
lineChartView1.TopAxisColor = [0.0, 0.0, 0.0]
lineChartView1.TopAxisTitleFontFamily = 'Arial'
lineChartView1.TopAxisTitleFontFile = ''
lineChartView1.TopAxisTitleFontSize = 18
lineChartView1.TopAxisTitleBold = 1
lineChartView1.TopAxisTitleItalic = 0
lineChartView1.TopAxisTitleColor = [0.0, 0.0, 0.0]
lineChartView1.TopAxisLogScale = 0
lineChartView1.TopAxisUseCustomRange = 0
lineChartView1.TopAxisRangeMinimum = 0.0
lineChartView1.TopAxisRangeMaximum = 1.0
lineChartView1.ShowTopAxisLabels = 1
lineChartView1.TopAxisLabelNotation = 'Mixed'
lineChartView1.TopAxisLabelPrecision = 2
lineChartView1.TopAxisUseCustomLabels = 0
lineChartView1.TopAxisLabels = []
lineChartView1.TopAxisLabelFontFamily = 'Arial'
lineChartView1.TopAxisLabelFontFile = ''
lineChartView1.TopAxisLabelFontSize = 12
lineChartView1.TopAxisLabelBold = 0
lineChartView1.TopAxisLabelItalic = 0
lineChartView1.TopAxisLabelColor = [0.0, 0.0, 0.0]

# show data in view
plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1)

# trace defaults for the display properties.
plotOverLine1Display_1.CompositeDataSetIndex = [0]
plotOverLine1Display_1.AttributeType = 'Point Data'
plotOverLine1Display_1.UseIndexForXAxis = 0
plotOverLine1Display_1.XArrayName = 'arc_length'
plotOverLine1Display_1.SeriesVisibility = ['Temperature']
plotOverLine1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'Temperature', 'Temperature', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'Temperature', '0.89', '0.1', '0.11', 'vtkValidPointMask', '0.22', '0.49', '0.72', 'Points_X', '0.3', '0.69', '0.29', 'Points_Y', '0.6', '0.31', '0.64', 'Points_Z', '1', '0.5', '0', 'Points_Magnitude', '0.65', '0.34', '0.16']
plotOverLine1Display_1.SeriesPlotCorner = ['arc_length', '0', 'Temperature', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display_1.SeriesLabelPrefix = ''
plotOverLine1Display_1.SeriesLineStyle = ['arc_length', '1', 'Temperature', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine1Display_1.SeriesLineThickness = ['arc_length', '2', 'Temperature', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'Temperature', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView1, layout=layout1, hint=0)

# update the view to ensure updated data information
lineChartView1.Update()

# save data
SaveData('/hpcwork/pn744180/parallel_computing/Homework_2/code/test/unrolled/T_finest.csv', proxy=plotOverLine1, WriteTimeSteps=0,
    Filenamesuffix='_%d',
    Precision=5,
    FieldDelimiter=',',
    UseScientificNotation=0,
    FieldAssociation='Points',
    AddMetaData=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 0.14142135623730953

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).