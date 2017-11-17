/*
 * HprfRfHammer2DGeneticSolver.java
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Comparable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Random;

import com.comsol.model.*;
import com.comsol.model.util.*;
/*
import org.jgap.Chromosome;
import org.jgap.FitnessFunction;
*/

/** Model exported on Jun 16 2015, 10:59 by COMSOL 4.4.0.195. */
// public class HprfRfHammer2DGeneticSolver extends FitnessFunction {
public class HprfRfHammer2DGeneticSolver {
  private Model model;
  private float[] times = null;
  // private ArrayList<float[]> data = new ArrayList<float[]>();
  private float[][][] data = null;
  private String tempFileDir;
  private String tempFileBase;
  private static final int PARAMETER_COUNT = 10;
  private static final int SWEEP_PARAMETER_COUNT = 9;
  private ArrayList<String> parameterSets = new ArrayList<String>();
  private static final String PARAM_STRING_TEMPLATE
    = "00000000000000000000000000000000000000000000000000000000000000000000000000";
  private static final String[] PARAM_STRING_TEMPLATES = new String[] {
    "000000", "000000", "000000", "000000000000", "00000000000",
    "0000", "000000000000", "00000000000", "00000", "00"};
  private static final String[] MIN_MAX_PARAM_STRINGS = new String[] {
    // "00000200000200000013400000000074000000000787009900000000044000000000889000",
    "070000070000000000134000000000740000000007870099000000000440000000001700000",
    "114300114300182880152000000000820000000008070125000000000490000000001900000"};
    /*
    "00000000000000000000001000000000000000000000001000000000000000000000000000010011010000000000000000000000000000000000000001110100000000000000000000000000000000000000011110000111000010011001000000000000000000000000000000000000010001000000000000000000000000000000000000001000100010010000",
    "00010001010000110000000000010001010000110000000000001001000101000100000000010101001000000000000000000000000000000000000010000010000000000000000000000000000000000000100000000111000000010010010100000000000000000000000000000000010010010000000000000000000000000000000000001000100110010000"};
    */
  private static final long[][] MIN_MAX_PARAMS = new long[][] {
    {070000, 070000, 0, 134000000000l, 74000000000l, 7870, 99000000000l, 44000000000l, 17000, 0},
    {114300, 114300, 182880, 152000000000l, 82000000000l, 8070, 125000000000l, 49000000000l, 19000, 0}};
    /*
    {2, 2, 0, 134000000000l, 74000000000l, 7870, 9900000000l, 44000000000l, 8890, 0},
    {114300, 114300, 182880, 152000000000l, 82000000000l, 8070, 12500000000l, 49000000000l, 8990, 30}};
    */
  private static final float SSR_MAX = 1.0e-3f;
  private int[] t0s = null;
  private int maxT0 = 0;
  private int duration = 0;

  private class Selection implements Comparable<Selection> {
    public Selection(final float fitness, final int index) {
      this.fitness = fitness;
      this.index = index;
    }
    public float fitness;
    public int index;
    public int compareTo(Selection obj) {
        float difference = obj.fitness - this.fitness;
        return Math.round(difference / Math.abs(difference));
    }
  }

  
  /****************** PUBLIC INTERFACE ******************/
  private HprfRfHammer2DGeneticSolver() {
    this.tempFileDir = "C:\\msys64\\tmp";
    this.tempFileBase = "temp_model";
  }
  
  public static HprfRfHammer2DGeneticSolver create() {
    HprfRfHammer2DGeneticSolver sim = new HprfRfHammer2DGeneticSolver();
    sim.createModel();
    sim.setDefaultParameters();
    sim.createFunctions();
    sim.createGeometry();
    sim.createMaterials();
    sim.createPhysics();
    sim.createMesh();
    sim.createStudy();
    sim.createSolver();
    sim.createBatch();
    sim.createPlots();
    
    return sim;
  }
  
  public static HprfRfHammer2DGeneticSolver load(File modelFile) throws java.io.IOException {
    HprfRfHammer2DGeneticSolver sim = new HprfRfHammer2DGeneticSolver();
    sim.model = ModelUtil.loadCopy(ModelUtil.uniquetag("Model"), modelFile.getAbsolutePath());
    
    return sim;
  }

  public File save() throws java.io.IOException {
    File tempFile = File.createTempFile(this.tempFileBase, ".mph", new File(this.tempFileDir));
    tempFile.deleteOnExit();
    this.model.save(tempFile.getAbsolutePath());

    return tempFile;
  }

  public void build() {
    buildGeometry();
    buildMaterials();
    buildMesh();
    buildPhysics();
    buildPlots();
  }

  public void run() {
    finalizeParameters();
    finalizeTimeSteps();
    
    System.out.println("Running simulation..");
    System.out.flush();
    // this.model.sol("sol1").runAll();
    this.model.batch("p1").run();
    
    String[] resultTags = this.model.result().tags();
    System.out.println("Results:");
    for (int index = 0; index < resultTags.length; ++index) {
      System.out.println("tag: " + resultTags[index]);
    }

    // The default solution selection is "all", so count the solutions
    ResultFeature plotGroup = this.model.result(resultTags[0]);
    final int SOLUTION_COUNT = parameterSets.size();
    System.out.println("# Solutions: " + SOLUTION_COUNT);

    // Get the number of plots we ought to have
    ResultFeature[] plots = new ResultFeature[] {this.model.result("pg3").feature("ptgr0"),
                                                 this.model.result("pg3").feature("ptgr5"),
                                                 this.model.result("pg3").feature("ptgr1")};
    // final int PLOT_COUNT = plot.selection().inputEntities().length;
    final int PLOT_COUNT = 3;

    // Extract the plot data
    this.data = null;
    boolean gotTimes = false;
    for (int solutionIndex = 0; solutionIndex < SOLUTION_COUNT; ++solutionIndex) {
      plotGroup.set("outersolnum", String.valueOf(solutionIndex+1));
      for (int plotIndex = 0; plotIndex < PLOT_COUNT; ++plotIndex) {
        int traceCount = 0;
        try {
          traceCount = plots[plotIndex].getGroups(0);  //Assume only one render group
        } catch (Exception e) {
          System.out.println(e.getMessage());
        }
        System.out.println("# Traces: " + traceCount);
        if (traceCount == 0) {
          continue;
        }
        if (!gotTimes) {
          this.times = plots[plotIndex].getVertices(0, 0)[0];
          gotTimes = true;

          this.data = new float[SOLUTION_COUNT][PLOT_COUNT][this.times.length];
          // initialize data to zeros
          for (int i = 0; i < SOLUTION_COUNT; ++i) {
            for (int j = 0; j < PLOT_COUNT; ++j) {
              for (int k = 0; k < this.times.length; ++k) {
                this.data[i][j][k] = 0;
              }
            }
          }
        }
        float[] samples = plots[plotIndex].getData(0, 0, "Height");
        System.arraycopy(normalizeSignal(samples), 0,
                         this.data[solutionIndex][plotIndex], 0,
                         samples.length);
        /*
        System.out.println("Solution " + (solutionIndex+1) + " Pressure Data (Pa):");
        for (int index=0; index < samples.length; ++index) {
          System.out.print(this.data[solutionIndex][plotIndex][index] + " ");
        }
        System.out.println();
        */
      }
    }
  }
/*
  public Model getModel() {
    return this.model;
  }
  
  public void setModel(Model model) {
    this.model = model;
  }
*/
  public float[][][] getData() {
    return this.data;
  }
  
  public float[] getTimes() {
    return this.times;
  }
  
  public void setDuration(int duration) {
    this.duration = duration;
  }

  public void clearParams() {
    this.parameterSets.clear();
  }

  public void addParamsFromString(String params) {
    this.parameterSets.add(params);
  }

  public String getParamsAsString(int paramSetIndex) {
    return this.parameterSets.get(paramSetIndex);
  }

  public float[] getParams(int paramSetIndex) {
    String paramsString = this.parameterSets.get(paramSetIndex);
    String[] genes = getGenes(paramsString);
    float[] parameters = new float[genes.length];
    for (int paramIndex = 0; paramIndex < parameters.length; ++paramIndex) {
      float parameter = 0.0f;
      if (paramIndex < 3) {
        parameter = Float.parseFloat("0." + genes[paramIndex]);
      } else {
        parameter = Float.parseFloat(genes[paramIndex]);
      }
      if (paramIndex == 2) {
        parameter -= MIN_MAX_PARAMS[1][2]/2 * 1e-6;  // subtract half the maximum for parameter 2
      }
      parameters[paramIndex] = parameter;
    }

    return parameters;
  }

  private static String binaryEncodeParamsString(String paramString) {
    StringBuffer encodedParamString = new StringBuffer(paramString.length() * 4);
    for (int decDigitIndex = 0; decDigitIndex < paramString.length(); ++decDigitIndex) {
      char decDigitValue = (char) (paramString.charAt(decDigitIndex) - '0');
      for (int binDigitIndex = 3; binDigitIndex >= 0; --binDigitIndex) {
        char bit = (char) ((decDigitValue & (1 << binDigitIndex)) >> binDigitIndex);
        encodedParamString.append((char)(bit + '0'));
      }
    }
    
    return encodedParamString.toString();
  }
  
  private static String binaryDecodeParamsString(String binaryString) {
    final int decDigitCount = binaryString.length() / 4;
    StringBuffer paramString = new StringBuffer(decDigitCount);
    for (int decDigitIndex = 0; decDigitIndex < decDigitCount; ++decDigitIndex) {
      int binDigitStartIndex = 4 * decDigitIndex;
      String decDigitBinString = binaryString.substring(binDigitStartIndex, binDigitStartIndex + 4);
      char decDigitValue = 0;
      for (int binDigitIndex = 0; binDigitIndex < 4; ++binDigitIndex) {
        char binDigit = (char) (decDigitBinString.charAt(binDigitIndex) - '0');
        decDigitValue = (char) (decDigitValue | (binDigit << (3-binDigitIndex)));
      }
      decDigitValue += '0';
      paramString.append(decDigitValue);
    }
    
    return paramString.toString();
  }

  /****************** PRIVATE METHODS ******************/
  private void createModel() {
    System.out.print("Creating model...");
    model = ModelUtil.create("Model");

    model.modelPath("C:\\Users\\plane\\Desktop\\COMSOL");

    model.name("HPRF RF Hammer 2D Export.mph");

    System.out.print("...");
    model.modelNode().create("comp1");
    System.out.println("Done");
  }

  private void createGeometry() {
    System.out.print("Creating geometry");
    System.out.flush();
    model.geom().create("geom1", 2);
    model.geom("geom1").axisymmetric(true);
    System.out.print(".");
    model.geom("geom1").feature().create("r12", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r4", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r1", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r2", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r5", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("dif1", "Difference");
    System.out.print(".");
    model.geom("geom1").feature().create("r3", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r6", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("dif2", "Difference");
    System.out.print(".");
    model.geom("geom1").feature().create("pol1", "Polygon");
    System.out.print(".");
    model.geom("geom1").feature().create("c1", "Circle");
    System.out.print(".");
    model.geom("geom1").feature().create("c2", "Circle");
    System.out.print(".");
    model.geom("geom1").feature().create("r7", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r10", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r14", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r17", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("r15", "Rectangle");
    System.out.print(".");
    model.geom("geom1").feature().create("b1", "BezierPolygon");
    System.out.print(".");
    model.geom("geom1").feature().create("pt1", "Point");
    System.out.print(".");
    model.geom("geom1").feature().create("pt2", "Point");
    System.out.println(".Done");
    model.geom("geom1").feature().create("pt3", "Point");
    model.geom("geom1").feature("r12").name("Environment");
    model.geom("geom1").feature("r12")
         .set("size", new String[]{"OuterRadius+EndHeight", "InnerHeight+4*EndHeight"});
    model.geom("geom1").feature("r12")
         .set("pos", new String[]{"0", "-InnerHeight/2-2*EndHeight"});
    model.geom("geom1").feature("r4").name("Cavity");
    model.geom("geom1").feature("r4")
         .set("size", new String[]{"15.24[cm]", "InnerHeight+2*EndHeight"});
    model.geom("geom1").feature("r4")
         .set("pos", new String[]{"0", "-InnerHeight/2-EndHeight"});
    model.geom("geom1").feature("r1").name("Resonance Chamber");
    model.geom("geom1").feature("r1")
         .set("size", new String[]{"InnerRadius", "InnerHeight"});
    model.geom("geom1").feature("r1")
         .set("pos", new String[]{"0", "-InnerHeight/2"});
    model.geom("geom1").feature("r2").name("Upstream Button Shim");
    model.geom("geom1").feature("r2")
         .set("size", new String[]{"ButtonRadius", "ShimHeight"});
    model.geom("geom1").feature("r2")
         .set("pos", new String[]{"0", "-InnerHeight/2"});
    model.geom("geom1").feature("r5").name("Upstream Button Shim Gap");
    model.geom("geom1").feature("r5")
         .set("size", new String[]{"ShimGapRadius", "ShimGapHeight"});
    model.geom("geom1").feature("r5")
         .set("pos", new String[]{"0", "-InnerHeight/2"});
    model.geom("geom1").feature("dif1").selection("input")
         .set(new String[]{"r2"});
    model.geom("geom1").feature("dif1").selection("input2")
         .set(new String[]{"r5"});
    model.geom("geom1").feature("r3").name("Downstream Button Shim");
    model.geom("geom1").feature("r3")
         .set("size", new String[]{"ButtonRadius", "ShimHeight"});
    model.geom("geom1").feature("r3")
         .set("pos", new String[]{"0", "InnerHeight/2-ShimHeight"});
    model.geom("geom1").feature("r6").name("Downstream Button Shim Gap");
    model.geom("geom1").feature("r6")
         .set("size", new String[]{"ShimGapRadius", "ShimGapHeight"});
    model.geom("geom1").feature("r6")
         .set("pos", new String[]{"0", "InnerHeight/2-ShimGapHeight"});
    model.geom("geom1").feature("dif2").selection("input")
         .set(new String[]{"r3"});
    model.geom("geom1").feature("dif2").selection("input2")
         .set(new String[]{"r6"});
    model.geom("geom1").feature("pol1")
         .set("table", new String[][]{{"ShimGapRadius+.15875[cm]", "-InnerHeight/2"}, {"ShimGapRadius", "-InnerHeight/2+ShimGapHeight"}});
    model.geom("geom1").feature("pol1").set("type", "open");
    model.geom("geom1").feature("pol1").set("source", "table");
    model.geom("geom1").feature("c1").name("Upstream Button");
    model.geom("geom1").feature("c1").set("angle", "90");
    model.geom("geom1").feature("c1").set("r", "ButtonRadius");
    model.geom("geom1").feature("c1")
         .set("pos", new String[]{"0", "-InnerHeight/2+ShimHeight"});
    model.geom("geom1").feature("c2").name("Downstream Button");
    model.geom("geom1").feature("c2").set("angle", "90");
    model.geom("geom1").feature("c2").set("r", "ButtonRadius");
    model.geom("geom1").feature("c2").set("rot", "270");
    model.geom("geom1").feature("c2")
         .set("pos", new String[]{"0", "InnerHeight/2-ShimHeight"});
    model.geom("geom1").feature("r7").name("Beam Hole");
    model.geom("geom1").feature("r7")
         .set("size", new String[]{"0.00508", "3.4925 [cm]"});
    model.geom("geom1").feature("r7")
         .set("pos", new String[]{"0", "-InnerHeight/2-EndHeight"});
    model.geom("geom1").feature("r10").name("Upstream Button Bolt Bore");
    model.geom("geom1").feature("r10")
         .set("size", new String[]{"0.00508", "4[cm]"});
    model.geom("geom1").feature("r10")
         .set("pos", new String[]{"0", "-InnerHeight/2-1.27[cm]"});
    model.geom("geom1").feature("r14").name("Upstream Button Bolt");
    model.geom("geom1").feature("r14")
         .set("size", new String[]{".4445 [cm]", "2.54 [cm]"});
    model.geom("geom1").feature("r14")
         .set("pos", new String[]{".508 [cm]", "-InnerHeight/2-1.27[cm]"});
    model.geom("geom1").feature("r17")
         .name("Upstream Button Bolt Drill Space");
    model.geom("geom1").feature("r17")
         .set("size", new String[]{".23335 [cm]", ".23335 [cm]"});
    model.geom("geom1").feature("r17")
         .set("pos", new String[]{".508 [cm]", "-InnerHeight/2+1.27 [cm]"});
    model.geom("geom1").feature("r15").name("Downstream Button Bolt");
    model.geom("geom1").feature("r15")
         .set("size", new String[]{".9525 [cm]", "2.54 [cm]"});
    model.geom("geom1").feature("r15")
         .set("pos", new String[]{"0", "InnerHeight/2-1.27[cm]"});
    model.geom("geom1").feature("b1").set("degree", new String[]{"1"});
    model.geom("geom1").feature("b1")
         .set("p", new String[][]{{"0.00508", "0"}, {"-0.01334", "-0.011"}});
    model.geom("geom1").feature("b1").set("w", new String[]{"1", "1"});
    model.geom("geom1").feature("pt1").name("S1");
    model.geom("geom1").feature("pt1").setIndex("p", "OuterRadius", 0, 0);
    model.geom("geom1").feature("pt1").setIndex("p", "z_S1", 1, 0);
    model.geom("geom1").feature("pt2").name("S0");
    model.geom("geom1").feature("pt2").setIndex("p", "r_S0", 0, 0);
    model.geom("geom1").feature("pt2")
         .setIndex("p", "-InnerHeight/2-EndHeight", 1, 0);
    model.geom("geom1").feature("pt3").name("S5");
    model.geom("geom1").feature("pt3").setIndex("p", "r_S5", 0, 0);
    model.geom("geom1").feature("pt3")
         .setIndex("p", "InnerHeight/2+EndHeight", 1, 0);
  }

  private void createFunctions() {
    System.out.print("Creating functions");
    model.func().create("pw1", "Piecewise");
    System.out.print(".");
    model.func().create("an1", "Analytic");
    System.out.print(".");
    model.func().create("an8", "Analytic");
    System.out.print(".");
    model.func().create("an9", "Analytic");
    System.out.print(".");
    model.func().create("an10", "Analytic");
    System.out.print(".");
    model.func("pw1")
         .set("pieces", new String[][]{{"0", "t_p/3", "log(1+3*(exp(1)-1)/t_p*t)"}, {"t_p/3", "3*t_p/4", "1"}, {"3*t_p/4", "t_p", "-log(1-4*(1-exp(-1)) + 4*(1-exp(-1))/t_p*t)"}});
    model.func("pw1").set("argunit", "s");
    model.func("pw1").set("funcname", "RFEnvelope");
    model.func("pw1").set("arg", "t");
    model.func("an1")
         .set("args", new String[]{"x", "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", 
         "a9", "a10", "a11", "a12"});
    model.func("an1")
         .set("expr", "a0+a1*x+a2*x^2+a3*x^3+a4*x^4+a5*x^5+a6*x^6+a7*x^7+a8*x^8+a9*x^9+a10*x^10+a11*x^11+a12*x^12");
    model.func("an1")
         .set("plotargs", new String[][]{{"x", "0", "1"}, {"a0", "0", "1"}, {"a1", "0", "1"}, {"a2", "0", "1"}, {"a3", "0", "1"}, {"a4", "0", "1"}, {"a5", "0", "1"}, {"a6", "0", "1"}, {"a7", "0", "1"}, {"a8", "0", "1"}, 
         {"a9", "0", "1"}, {"a10", "0", "1"}, {"a11", "0", "1"}, {"a12", "0", "1"}});
    model.func("an1").set("funcname", "polynomial");
    model.func("an8").set("args", new String[]{"rho"});
    model.func("an8").set("argunit", "cm");
    model.func("an8")
         .set("expr", "polynomial(rho,-240.178811905,-278.60505448,-178.332730546,72.9329677377,174.862318937,-146.980411248,53.7093343214,-11.4886388004,1.56357151204,-0.137733463426,0.00763183423149,-0.000242374455041,3.36971433768e-06)");
    model.func("an8")
         .set("plotargs", new String[][]{{"rho", "ButtonRadius", "InnerRadius"}});
    model.func("an8").set("funcname", "f_wall");
    model.func("an8").set("fununit", "N/m^2");
    model.func("an9").set("args", new String[]{"z"});
    model.func("an9").set("argunit", "cm");
    model.func("an9")
         .set("expr", "polynomial(z,-84.0445115747,9.21244614012e-05,-0.0257325876635,-0.00010904929815,0.00109329570891,4.34025393603e-05,2.64423467733e-05,-6.6999234955e-06,-6.44621761842e-06,4.29603467484e-07,3.99598700887e-07,-9.44057804598e-09,-8.75366850967e-09)");
    model.func("an9")
         .set("plotargs", new String[][]{{"z", "-InnerHeight/2", "InnerHeight/2"}});
    model.func("an9").set("funcname", "f_circ");
    model.func("an9").set("fununit", "N/m^2");
    model.func("an10").set("args", new String[]{"z"});
    model.func("an10").set("argunit", "cm");
    model.func("an10")
         .set("expr", "polynomial(z,15758.5245901,-19.6518218233,-9618.87504383,15.8366516429,2607.34389555,-4.31630431287,-372.886909702,0.521172056158,28.8134963162,-0.028761044486,-1.13644437236,0.000592426367832,0.0179375462308)");
    model.func("an10")
         .set("plotargs", new String[][]{{"z", "-InnerHeight/2", "InnerHeight/2"}});
    model.func("an10").set("funcname", "f_button");
    model.func("an10").set("fununit", "N/m^2");
    System.out.println("Done");
  }

  private void createMaterials() {
    System.out.print("Creating materials");
    Material nitrogen = model.material().create("mat4");
    System.out.print(".");
    nitrogen.propertyGroup("def").func().create("eta", "Piecewise");
    System.out.print(".");
    nitrogen.propertyGroup("def").func().create("Cp", "Piecewise");
    System.out.print(".");
    nitrogen.propertyGroup("def").func().create("rho", "Analytic");
    System.out.print(".");
    nitrogen.propertyGroup("def").func().create("k", "Piecewise");
    System.out.print(".");

    Material copper = model.material().create("mat2");
    copper.propertyGroup().create("Enu", "Young's modulus and Poisson's ratio");
    copper.propertyGroup().create("linzRes", "Linearized resistivity");
    copper.propertyGroup()
         .create("CpCs", "Pressure-wave and shear-wave speeds");
    copper.propertyGroup()
         .create("KG", "Bulk modulus and shear modulus");
    System.out.print(".");

    Material steel = model.material().create("mat3");
    steel.propertyGroup().create("comp1", "User-defined property group");
    steel.propertyGroup().create("CpCs", "Pressure-wave and shear-wave speeds");
    steel.propertyGroup().create("KG", "Bulk modulus and shear modulus");
    System.out.print(".");

    Material air = model.material().create("mat5");
    air.propertyGroup("def").func().create("eta", "Piecewise");
    air.propertyGroup("def").func().create("Cp", "Piecewise");
    air.propertyGroup("def").func().create("rho", "Analytic");
    air.propertyGroup("def").func().create("k", "Piecewise");
    air.propertyGroup("def").func().create("cs", "Analytic");
    System.out.print(".");

    nitrogen.name("Nitrogen");
    nitrogen.propertyGroup("def").func("eta")
      .set("pieces", new String[][]{{"200.0", "1200.0", "1.77230303E-6+6.27427545E-8*T^1-3.47278555E-11*T^2+1.01243201E-14*T^3"}});
    nitrogen.propertyGroup("def").func("eta").set("arg", "T");
    nitrogen.propertyGroup("def").func("Cp")
      .set("pieces", new String[][]{{"200.0", "1200.0", "1088.22121-0.365941919*T^1+7.88715035E-4*T^2-3.749223E-7*T^3+3.17599068E-11*T^4"}});
    nitrogen.propertyGroup("def").func("Cp").set("arg", "T");
    nitrogen.propertyGroup("def").func("rho").set("args", new String[]{"pA", "T"});
    nitrogen.propertyGroup("def").func("rho").set("expr", "pA*0.02801/8.314/T");
    nitrogen.propertyGroup("def").func("rho").set("dermethod", "manual");
    nitrogen.propertyGroup("def").func("rho")
      .set("plotargs", new String[][]{{"pA", "0", "1"}, {"T", "0", "1"}});
    nitrogen.propertyGroup("def").func("rho")
         .set("argders", new String[][]{{"pA", "d(pA*0.02801/8.314/T,pA)"}, {"T", "d(pA*0.02801/8.314/T,T)"}});
    nitrogen.propertyGroup("def").func("k")
         .set("pieces", new String[][]{{"200.0", "1200.0", "3.6969697E-4+9.74353924E-5*T^1-4.07587413E-8*T^2+7.68453768E-12*T^3"}});
    nitrogen.propertyGroup("def").func("k").set("arg", "T");
    nitrogen.propertyGroup("def")
         .set("dynamicviscosity", "eta(T[1/K])[Pa*s]");
    nitrogen.propertyGroup("def")
         .set("ratioofspecificheat", "1.4");
    nitrogen.propertyGroup("def")
         .set("heatcapacity", "Cp(T[1/K])[J/(kg*K)]");
    nitrogen.propertyGroup("def").set("density", "37.5");
    nitrogen.propertyGroup("def")
         .set("thermalconductivity", new String[]{"k(T[1/K])[W/(m*K)]", "0", "0", "0", "k(T[1/K])[W/(m*K)]", "0", "0", "0", "k(T[1/K])[W/(m*K)]"});
    nitrogen.propertyGroup("def").set("soundspeed", "c_N2");
    nitrogen.propertyGroup("def").addInput("temperature");
    nitrogen.propertyGroup("def").addInput("pressure");
    System.out.print(".");

    copper.name("Copper");
    copper.propertyGroup("def")
         .set("relpermeability", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    copper.propertyGroup("def")
         .set("electricconductivity", new String[]{"5.998e7[S/m]", "0", "0", "0", "5.998e7[S/m]", "0", "0", "0", "5.998e7[S/m]"});
    copper.propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"17e-6[1/K]", "0", "0", "0", "17e-6[1/K]", "0", "0", "0", "17e-6[1/K]"});
    copper.propertyGroup("def")
         .set("heatcapacity", "385[J/(kg*K)]");
    copper.propertyGroup("def")
         .set("relpermittivity", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    copper.propertyGroup("def").set("density", "rho_Cu");
    copper.propertyGroup("def")
         .set("thermalconductivity", new String[]{"400[W/(m*K)]", "0", "0", "0", "400[W/(m*K)]", "0", "0", "0", "400[W/(m*K)]"});
    copper.propertyGroup("def").set("soundspeed", "4760");
    copper.propertyGroup("Enu").set("youngsmodulus", "117e9");
    copper.propertyGroup("Enu").set("poissonsratio", "0.35");
    copper.propertyGroup("linzRes").set("rho0", "");
    copper.propertyGroup("linzRes").set("alpha", "");
    copper.propertyGroup("linzRes").set("Tref", "");
    copper.propertyGroup("linzRes")
         .set("rho0", "1.72e-8[ohm*m]");
    copper.propertyGroup("linzRes")
         .set("alpha", "0.0039[1/K]");
    copper.propertyGroup("linzRes").set("Tref", "298[K]");
    copper.propertyGroup("linzRes").addInput("temperature");
    copper.propertyGroup("CpCs").set("cp", "");
    copper.propertyGroup("CpCs").set("cs", "");
    copper.propertyGroup("CpCs").set("cp", "c_p_Cu");
    copper.propertyGroup("CpCs").set("cs", "c_s_Cu");
    copper.propertyGroup("KG").set("K", "");
    copper.propertyGroup("KG").set("G", "");
    copper.propertyGroup("KG").set("K", "K_Cu");
    copper.propertyGroup("KG").set("G", "G_Cu");
    System.out.print(".");

    steel.name("316 Stainless Steel");
    steel.propertyGroup("def").set("density", "rho_SS");
    steel.propertyGroup("def").set("soundspeed", "5790");
    steel.propertyGroup("def")
         .set("relpermeability", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    steel.propertyGroup("def")
         .set("electricconductivity", new String[]{"1.42e6", "0", "0", "0", "1.42e6", "0", "0", "0", "1.42e6"});
    steel.propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"1.6e-5", "0", "0", "0", "1.6e-5", "0", "0", "0", "1.6e-5"});
    steel.propertyGroup("def").set("heatcapacity", "500");
    steel.propertyGroup("def")
         .set("relpermittivity", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    steel.propertyGroup("def")
         .set("thermalconductivity", new String[]{"16.2", "0", "0", "0", "16.2", "0", "0", "0", "16.2"});
    steel.propertyGroup("def")
         .set("youngsmodulus", "197.5e9");
    steel.propertyGroup("def").set("poissonsratio", "0.30");
    steel.propertyGroup("def").set("lossfactor", ".04");
    steel.propertyGroup("CpCs").set("cp", "");
    steel.propertyGroup("CpCs").set("cs", "");
    steel.propertyGroup("CpCs").set("cp", "c_p_SS");
    steel.propertyGroup("CpCs").set("cs", "c_s_SS");
    steel.propertyGroup("KG").set("K", "");
    steel.propertyGroup("KG").set("G", "");
    steel.propertyGroup("KG").set("K", "K_SS");
    steel.propertyGroup("KG").set("G", "G_SS");
    System.out.print(".");

    air.name("Air");
    air.propertyGroup("def").func("eta")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4"}});
    air.propertyGroup("def").func("eta").set("arg", "T");
    air.propertyGroup("def").func("Cp")
         .set("pieces", new String[][]{{"200.0", "1600.0", "1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4"}});
    air.propertyGroup("def").func("Cp").set("arg", "T");
    air.propertyGroup("def").func("rho")
         .set("args", new String[]{"pA", "T"});
    air.propertyGroup("def").func("rho")
         .set("expr", "pA*0.02897/8.314/T");
    air.propertyGroup("def").func("rho")
         .set("dermethod", "manual");
    air.propertyGroup("def").func("rho")
         .set("plotargs", new String[][]{{"pA", "0", "1"}, {"T", "0", "1"}});
    air.propertyGroup("def").func("rho")
         .set("argders", new String[][]{{"pA", "d(pA*0.02897/8.314/T,pA)"}, {"T", "d(pA*0.02897/8.314/T,T)"}});
    air.propertyGroup("def").func("k")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4"}});
    air.propertyGroup("def").func("k").set("arg", "T");
    air.propertyGroup("def").func("cs")
         .set("args", new String[]{"T"});
    air.propertyGroup("def").func("cs")
         .set("expr", "sqrt(1.4*287*T)");
    air.propertyGroup("def").func("cs")
         .set("dermethod", "manual");
    air.propertyGroup("def").func("cs")
         .set("plotargs", new String[][]{{"T", "0", "1"}});
    air.propertyGroup("def").func("cs")
         .set("argders", new String[][]{{"T", "d(sqrt(1.4*287*T),T)"}});
    air.propertyGroup("def")
         .set("dynamicviscosity", "eta(T[1/K])[Pa*s]");
    air.propertyGroup("def")
         .set("ratioofspecificheat", "1.4");
    air.propertyGroup("def")
         .set("electricconductivity", new String[]{"0[S/m]", "0", "0", "0", "0[S/m]", "0", "0", "0", "0[S/m]"});
    air.propertyGroup("def")
         .set("heatcapacity", "Cp(T[1/K])[J/(kg*K)]");
    air.propertyGroup("def")
         .set("density", "rho(pA[1/Pa],T[1/K])[kg/m^3]");
    air.propertyGroup("def")
         .set("thermalconductivity", new String[]{"k(T[1/K])[W/(m*K)]", "0", "0", "0", "k(T[1/K])[W/(m*K)]", "0", "0", "0", "k(T[1/K])[W/(m*K)]"});
    air.propertyGroup("def")
         .set("soundspeed", "cs(T[1/K])[m/s]");
    air.propertyGroup("def").addInput("temperature");
    air.propertyGroup("def").addInput("pressure");
    System.out.println("Done");
  }

  private void createPhysics() {
    System.out.print("Creating physics");
    model.physics().create("astd", "TransientAcousticSolid", "geom1");
    System.out.print(".");
    model.physics("astd").field("pressure").field("p2");
    model.physics("astd").feature("dlemm1").feature()
         .create("dmp1", "Damping", 2);
    model.physics("astd").feature().create("dlemm2", "LinearElasticModel", 2);
    System.out.print(".");
    model.physics("astd").feature("dlemm2").feature()
         .create("dmp1", "Damping", 2);
    model.physics("astd").feature().create("pwr1", "PlaneWaveRadiation", 1);
    System.out.print(".");
    model.physics("astd").feature().create("bndl1", "BoundaryLoad", 1);
    System.out.print(".");
    model.physics("astd").feature().create("bndl2", "BoundaryLoad", 1);
    System.out.print(".");
    model.physics("astd").feature().create("bndl3", "BoundaryLoad", 1);
    System.out.print(".");

    model.physics("astd").feature("tpam1").set("editModelInputs", "1");
    model.physics("astd").feature("tpam1").set("minput_pressure", "1 [atm]");
    model.physics("astd").feature("dlemm1").set("IsotropicOption", "KG");
    model.physics("astd").feature("dlemm1").feature("dmp1")
         .set("beta_dK", "0");
         // .set("beta_dK", "1e-8");
    model.physics("astd").feature("dlemm2").set("IsotropicOption", "KG");
    model.physics("astd").feature("dlemm2").feature("dmp1")
         .set("beta_dK", "0");
         // .set("beta_dK", "5e-4");
    model.physics("astd").feature("bndl1")
         .set("FperArea", new String[][]{{"0"}, {"0"}, {"-z/abs(z)*RFEnvelope(t)*f_wall(r)"}});
    model.physics("astd").feature("bndl1").set("coordinateSystem", "sys1");
    model.physics("astd").feature("bndl1")
         .name("End Flange RF Hammer Envelope Load");
    model.physics("astd").feature("bndl2")
         .set("FperArea", new String[][]{{"0"}, {"0"}, {"-RFEnvelope(t)*f_circ(z)"}});
    model.physics("astd").feature("bndl2").set("coordinateSystem", "sys1");
    model.physics("astd").feature("bndl2")
         .name("Center Flange RF Hammer Envelope Load");
    model.physics("astd").feature("bndl3")
         .set("FperArea", new String[][]{{"0"}, {"0"}, {"RFEnvelope(t)*f_button(z)"}});
    model.physics("astd").feature("bndl3").set("coordinateSystem", "sys1");
    model.physics("astd").feature("bndl3")
         .name("Button RF Hammer Envelope Load");
    System.out.println("Done");
  }

  private void createMesh() {
    System.out.print("Creating mesh");
    model.mesh().create("mesh1", "geom1");
    System.out.print(".");
    model.mesh("mesh1").feature().create("ftri1", "FreeTri");
    System.out.print(".");
    model.mesh("mesh1").feature("size").set("hauto", 3);
    System.out.println("Done");
  }
  
  private void createStudy() {
    System.out.print("Creating study");
    Study study = this.model.study().create("std1");
    System.out.print(".");
    study.feature().create("time", "Transient");
    System.out.print(".");

    study.feature("time").set("initstudyhide", "on");
    study.feature("time").set("initsolhide", "on");
    study.feature("time").set("notstudyhide", "on");
    study.feature("time").set("notsolhide", "on");

    //study.feature("time").set("tlist", "range(0,2[us],500[us])");
    study.feature("time").set("tlist", "range(0,2[us],500[us])");
    System.out.print(".");

    study.feature().create("param", "Parametric");
    study.feature("param")
         .set("plistarr", new String[]{"0.097967", "0.100000", "0.000000", "134000000000", "75000000000", "8070", "99000000000", "44000000000", "8990"});
    study.feature("param")
         .set("pname", new String[]{"r_S0", "r_S5", "z_S1", "K_SS", "G_SS", "rho_SS", "K_Cu", "G_Cu", "rho_Cu"});
    System.out.println("Done");
  }

  private void createSolver() {
    System.out.print("Creating solver");
    this.model.sol().create("sol1");
    System.out.print(".");
    this.model.sol("sol1").study("std1");
    this.model.sol("sol1").attach("std1");
    System.out.print(".");
    this.model.sol("sol1").feature().create("st1", "StudyStep");
    System.out.print(".");
    this.model.sol("sol1").feature().create("v1", "Variables");
    System.out.print(".");
    this.model.sol("sol1").feature().create("t1", "Time");
    System.out.print(".");
    this.model.sol("sol1").feature("t1").feature().create("fc1", "FullyCoupled");
    this.model.sol("sol1").feature("t1").feature().remove("fcDef");
    System.out.print(".");

    this.model.sol("sol1").feature("st1")
         .name("Compile Equations: Time Dependent");
    this.model.sol("sol1").feature("st1").set("studystep", "time");
    this.model.sol("sol1").feature("v1").set("control", "time");
    // this.model.sol("sol1").feature("t1").set("tlist", "range(0,2[us],500[us])");
    this.model.sol("sol1").feature("t1").set("tlist", "range(0,2[us],500[us])");
    this.model.sol("sol1").feature("t1").set("timemethod", "genalpha");
    this.model.sol("sol1").feature("t1").set("timestepgenalpha", "2e-6");
    this.model.sol("sol1").feature("t1").set("tstepsgenalpha", "manual");
    this.model.sol("sol1").feature("t1").set("control", "time");
    this.model.sol("sol1").feature("t1").set("rhoinf", "0.2");
    System.out.println("Done");

    this.model.sol().create("sol2");
    this.model.sol("sol2").study("std1");
    this.model.sol("sol2").feature().create("su1", "StoreSolution");
  }

  private void createBatch() {
    System.out.print("Creating batch");
    this.model.batch().create("p1", "Parametric");
    System.out.print(".");
    this.model.batch("p1").feature().create("so1", "Solutionseq");
    System.out.print(".");
    this.model.batch("p1").study("std1");
    this.model.batch("p1")
         .set("plistarr", new String[]{"0.097967", "0.100000", "0.000000", "134000000000", "75000000000", "8070", "99000000000", "44000000000", "8990"});
    this.model.batch("p1").set("err", true);
    this.model.batch("p1")
         .set("pname", new String[]{"r_S0", "r_S5", "z_S1", "K_SS", "G_SS", "rho_SS", "K_Cu", "G_Cu", "rho_Cu"});
    System.out.print(".");
    this.model.batch("p1").set("control", "param");
    this.model.batch("p1").attach("std1");
    this.model.batch("p1").feature("so1").set("psol", "sol2");
    this.model.batch("p1").feature("so1").set("seq", "sol1");
    System.out.println("Done");
  }

  private void createPlots() {
    System.out.print("Creating plots");
    this.model.result().create("pg3", "PlotGroup1D");
    System.out.print(".");
    this.model.result("pg3").name("1D Solid Plots");
    this.model.result("pg3").set("ylabel", "Pressure (N/m<sup>2</sup>)");
    this.model.result("pg3").set("manualgrid", "on");
    this.model.result("pg3").set("xlabel", "t (s)");
    this.model.result("pg3").set("xspacing", "2.0E-6");
    this.model.result("pg3").set("yspacing", "1.0E-12");
    this.model.result("pg3").set("xlabelactive", false);
    this.model.result("pg3").set("ylabelactive", false);
    System.out.print(".");
    this.model.result("pg3").set("data", "dset2");

    this.model.result("pg3").feature().create("ptgr0", "PointGraph");
    this.model.result("pg3").feature("ptgr0").set("xdatasolnumtype", "inner");
    ResultFeature ptgr0 = model.result("pg3").feature("ptgr0");
    ptgr0.name("S0 Acceleration");
    ptgr0.set("expr", "astd.u_ttZ");
    ptgr0.set("unit", "N/m^2");

    this.model.result("pg3").feature().create("ptgr5", "PointGraph");
    this.model.result("pg3").feature("ptgr5").set("xdatasolnumtype", "inner");
    ResultFeature ptgr5 = model.result("pg3").feature("ptgr5");
    ptgr5.name("S5 Acceleration");
    ptgr5.set("expr", "-astd.u_ttZ");
    ptgr5.set("unit", "N/m^2");
    
    this.model.result("pg3").feature().create("ptgr1", "PointGraph");
    this.model.result("pg3").feature("ptgr1").set("xdatasolnumtype", "inner");
    ResultFeature ptgr1 = model.result("pg3").feature("ptgr1");
    ptgr1.name("S1 Acceleration");
    ptgr1.set("expr", "-astd.u_ttR");
    ptgr1.set("unit", "N/m^2");
	/*
    model.result("pg3").feature().create("ptgr8", "PointGraph");
    model.result("pg3").feature("ptgr8").active(false);
    model.result("pg3").feature("ptgr8").name("S0 Z Acceleration");
    model.result("pg3").feature("ptgr8").set("expr", "3e-10*astd.u_ttZ");
    model.result("pg3").feature("ptgr8").set("unit", "m/s^2");
    model.result("pg3").feature("ptgr8").set("descr", "3e-10*astd.u_ttZ");
    model.result("pg3").feature("ptgr8").selection().set(new int[]{49});

    model.result("pg3").feature().create("ptgr9", "PointGraph");
    model.result("pg3").feature("ptgr9").active(false);
    model.result("pg3").feature("ptgr9").name("S5 Z Acceleration");
    model.result("pg3").feature("ptgr9").set("expr", "3e-10*astd.u_ttZ");
    model.result("pg3").feature("ptgr9").set("unit", "m/s^2");
    model.result("pg3").feature("ptgr9").set("descr", "3e-10*astd.u_ttZ");
    model.result("pg3").feature("ptgr9").selection().set(new int[]{50});

    model.result("pg3").feature().create("ptgr5", "PointGraph");
    model.result("pg3").feature("ptgr5").active(false);
    model.result("pg3").feature("ptgr5").name("S1 Pressure");
    model.result("pg3").feature("ptgr5").set("expr", "astd.pm");
    model.result("pg3").feature("ptgr5").set("unit", "N/m^2");
    model.result("pg3").feature("ptgr5").selection().set(new int[]{54});
	*/
    System.out.println("Done");
  }

  private void buildGeometry() {
    System.out.print("Building geometry...");
    model.geom("geom1").run();
    System.out.println("Done");
  }
  
  private void buildMaterials() {
    System.out.print("Building materials...");
    // FIXME: don't rely on entity numbers
    this.model.material("mat4").selection().set(new int[]{10, 14, 21, 24, 25});
    this.model.material("mat2").selection().set(new int[]{9, 11, 17, 18, 19, 22, 23});
    this.model.material("mat3").selection().set(new int[]{1, 2, 3, 12, 13, 14, 15, 16, 17, 18, 19});
    this.model.material("mat5").selection().set(new int[]{1, 2, 4, 5, 6, 7, 8, 20});
    System.out.println("Done");
  }

  private void buildMesh() {
    System.out.print("Building mesh...");
    this.model.mesh("mesh1").run();
    System.out.println("Done");
  }

  private void buildPhysics() {
    System.out.print("Building physics...");
    // Linear elastic domains 1
    model.physics("astd").feature("dlemm1").selection()
         .set(new int[]{3, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23});
    // Linear elastic domains 2
    model.physics("astd").feature("dlemm2").selection()
         .set(new int[]{9, 11, 22, 23});
    // Plane wave radiation boundaries
    model.physics("astd").feature("pwr1").selection()
         .set(new int[]{2, 33, 79});
    // End Flange RF Hammer Envelope Load boundaries
    model.physics("astd").feature("bndl1").selection().set(new int[]{71, 73});
    // Center Flange RF Hammer Envelope Load boundaries
    model.physics("astd").feature("bndl2").selection().set(new int[]{76});
    // Button RF Hammer Envelope Load boundaries
    model.physics("astd").feature("bndl3").selection()
         .set(new int[]{70, 72, 80, 81});
    System.out.println("Done");
  }

  private void buildPlots() {
    System.out.print("Building plots...");
    /**/
    SelectionFeature S0_sel = model.selection().create("S0", "Ball");
    S0_sel.set("entitydim", "0");
    S0_sel.set("posx", "r_S0");
    S0_sel.set("posy", "-InnerHeight/2-EndHeight");
    S0_sel.set("r","1[mm]");
    SelectionFeature S1_sel = model.selection().create("S1", "Ball");
    S1_sel.set("entitydim", "0");
    S1_sel.set("posx", "OuterRadius");
    S1_sel.set("posy", "z_S1");
    S1_sel.set("r","1[mm]");
    SelectionFeature S5_sel = model.selection().create("S5", "Ball");
    S5_sel.set("entitydim", "0");
    S5_sel.set("posx", "r_S5");
    S5_sel.set("posy", "InnerHeight/2+EndHeight");
    S5_sel.set("r","1[mm]");
    /*
    SelectionFeature microphones_sel = model.selection().create("microphones", "Union");
    model.selection("microphones").set("entitydim", "0");
    model.selection("microphones").set("input", new String[]{"S0", "S1", "S5"});
    */
    /*
    final int UPSTREAM_MIC = 50;
    final int DOWNSTREAM_MIC = 49;
    final int CIRCUMFERENCE_MIC = 54;
    */
    
    ResultFeature ptgr0 = model.result("pg3").feature("ptgr0");
    // ptgr0.selection().set(new int[]{UPSTREAM_MIC});
    ptgr0.selection().named("S0");
    
    ResultFeature ptgr5 = model.result("pg3").feature("ptgr5");
    // ptgr5.selection().set(new int[]{DOWNSTREAM_MIC});
    ptgr5.selection().named("S5");

    ResultFeature ptgr1 = model.result("pg3").feature("ptgr1");
    // ptgr1.selection().set(new int[]{CIRCUMFERENCE_MIC});
    ptgr1.selection().named("S1");
	/*
    model.result("pg3").feature().create("ptgr10", "PointGraph");
    model.result("pg3").feature("ptgr10").name("S5 Pressure");
    model.result("pg3").feature("ptgr10").set("expr", "astd.pm");
    model.result("pg3").feature("ptgr10").set("unit", "N/m^2");
    model.result("pg3").feature("ptgr10").selection().set(new int[]{49, 50});
	*/

	/*
    model.result().export().create("plot1", "Plot");
    model.result().export("plot1")
         .set("filename", "C:\\Users\\plane\\Desktop\\COMSOL\\Analysis\\RF Hammer\\I_Got_Gas.csv");
	*/
    System.out.println("Done");
  }

  private void setDefaultParameters() {
    model.param().set("InnerRadius", "11.43 [cm]");
    model.param().set("OuterRadius", "15.24 [cm]");
    model.param().set("InnerHeight", "8.128 [cm]");
    model.param().set("EndHeight", "5.08 [cm]");
    model.param().set("ButtonRadius", "2.54 [cm]");
    model.param().set("ShimHeight", "0.635 [cm]");
    model.param().set("ShimGapRadius", "2.2225 [cm]");
    model.param().set("ShimGapHeight", "0.15748 [cm]");
    model.param().set("c_p_SS", "5800 [m/s]", "5790 [m/s]");
    model.param().set("c_s_SS", "0.52 * c_p_SS", "3100 [m/s]");
    model.param().set("c_p_Cu", "4480 [m/s]", "4760 [m/s]");
    model.param().set("c_s_Cu", "0.52 * c_p_Cu", "2325 [m/s]");
    model.param().set("c_N2", "540 [m/s]");
    model.param().set("t_p", "30 [us]");
    model.param().set("bolt_disp", "0 [cm]");
    model.param().set("notch_height", ".6731");
    model.param().set("RFrac", "0.7");
    model.param().set("r_S0", "0.097967");  // "0.8571*InnerRadius"
    model.param().set("r_S5", "0.100000"); // 10 [cm]"
    model.param().set("z_S1", "0.000000");
    model.param().set("K_SS", "134000000000");  // Pa
    model.param().set("G_SS", "75000000000");  // Pa
    model.param().set("rho_SS", "8070");  // kg/m^3
    model.param().set("K_Cu", "99000000000");  // Pa
    model.param().set("G_Cu", "44000000000");  // Pa
    model.param().set("rho_Cu", "8990");  // kg/m^3
  }

  private void finalizeParameters() {
    StringBuffer[] paramValues = new StringBuffer[this.PARAMETER_COUNT];
    boolean initialized = false;
    Iterator<String> paramSet = this.parameterSets.iterator();
    this.t0s = new int[this.parameterSets.size()];
    int individualIndex = 0;
    this.maxT0 = 0;
    while (paramSet.hasNext()) {
      String paramsString = paramSet.next();
      String[] paramStrings = getGenes(paramsString);
      paramStrings[0] = "0." + paramStrings[0];  // r_S0
      paramStrings[1] = "0." + paramStrings[1];  // r_S5

      final double z_S1 = (Integer.parseInt(paramStrings[2]) - (MIN_MAX_PARAMS[1][2]/2)) * 1e-6;
      paramStrings[2] = String.valueOf(z_S1); // z_S1

      // Convert sweep parameter sets to CSV strings per parameter
      for (int index = 0; index < HprfRfHammer2DGeneticSolver.SWEEP_PARAMETER_COUNT; ++index) {
        if (!initialized) {
          paramValues[index] = new StringBuffer();
        } else {
          paramValues[index].append(",");
        }
        paramValues[index].append(paramStrings[index]);
      }
      initialized = true;

      // set t0 values for each individual, and keep track of the largest for setting sim duration
      this.t0s[individualIndex] = Integer.parseInt(paramsString.substring(72,74));  // time offset
      if (this.t0s[individualIndex] > this.maxT0) {
        this.maxT0 = this.t0s[individualIndex];
      }
      ++individualIndex;
    }

    // Set the CSV strings for each parameter
    String[] parameters = new String[this.PARAMETER_COUNT-1];
    for (int index = 0; index < parameters.length; ++index) {
      parameters[index] = paramValues[index].toString();
      System.out.print("\"" + parameters[index] + "\" ");
    }
    System.out.println();
    this.model.study("std1").feature("param").set("plistarr", parameters);
    this.model.batch("p1").set("plistarr", parameters);
  }
  
  private void finalizeTimeSteps() {
    String range = new String("range(0,2[us]," + (this.duration + this.maxT0) + "[us])");
    this.model.study("std1").feature("time").set("tlist", range);
    this.model.sol("sol1").feature("t1").set("tlist", range);
  }

  public String select(final float sumOfFitnesses, final float[] fitnesses) {
    Random rand = new Random();
    // System.out.println("Sum of Fitnesses: " + sumOfFitnesses);
    float target = Math.abs(sumOfFitnesses * rand.nextFloat());
    float partialSum = 0.0f;
    int index = 0;
    do {
      partialSum += Math.abs(fitnesses[index]);
      ++index;
    } while ((partialSum < target) && (index < fitnesses.length));
    
    return this.parameterSets.get(index-1);
  }

  /* Select how many copies of an individual will be used in the mating process via the
   * Stochastic Universal Sampling method.
   */
  public ArrayList<String> selectFittest(final float[] fitnesses) {
    float sumOfFitnesses = 0.0f;
    for (int index = 0; index < fitnesses.length; ++index) {
      sumOfFitnesses += fitnesses[index];
    }
    final float avgFitness = sumOfFitnesses / fitnesses.length;

    float partialSum = 0;
    Random rand = new Random();
    float fitnessPointer = rand.nextFloat();
    ArrayList<Selection> selections = new ArrayList<Selection>(fitnesses.length);
    for (int index = 0; index < fitnesses.length; ++index) {
      float reproductionRate = 0.0f;
      if (avgFitness > 0) {
        reproductionRate = fitnesses[index] / avgFitness;  // proportional selection only
      } else {
        reproductionRate = 1.0f;
      }
      System.out.println("Reproduction Rate: " + reproductionRate);
      partialSum += reproductionRate;
      while (partialSum > fitnessPointer) {
        // rawSelections.add(this.parameterSets.get(index));
        selections.add(new Selection(fitnesses[index], index));
        System.out.println("Initial Selected Param String: " + this.parameterSets.get(index));
        fitnessPointer += 1.0f;
      }
    }

    Collections.sort(selections);

    // Trade two low-fitness individuals for extra copies of the top two individuals (elite)
    int[] fittestIndicies = new int[] {selections.get(0).index, 0};
    for (Selection selection: selections) {
      if (selection.index != fittestIndicies[0]) {
        fittestIndicies[1] = selection.index;
        break;
      }
    }
    for (int index = 1; index >= 0; --index) {
      selections.remove(selections.size() - 1);
      selections.add(0, new Selection(fitnesses[fittestIndicies[index]], fittestIndicies[index]));
    }

    // Occasionally the number of selections will be off by 1. Fix it.
    if (selections.size() > fitnesses.length) {
      // drop one of the entries with the lowest reproduction rate
      selections.remove(fitnesses.length-1);
    } else if (selections.size() < fitnesses.length) {
      // add an extra copy of one of the entries with the highest reproduction rate
      selections.add(0, selections.get(0));
    }

    System.out.println("New size of selections: " + selections.size());

    // Randomly sample without replacement the sorted list of selections
    ArrayList<String> selectedParams = new ArrayList<String>(fitnesses.length);
    for (int index = 0; index < fitnesses.length; ++index) {
      final int lastIndex = fitnesses.length - index - 1;
      final int selectionIndex = Math.round(lastIndex * rand.nextFloat());
      final Selection selection = selections.remove(selectionIndex);
      String selectedParamString = this.parameterSets.get(selection.index);
      System.out.println("Selected Param String: " + selectedParamString);
      selectedParams.add(selectedParamString);
    }

    return selectedParams;
  }

  private static float[] scalingCoefficients(final float minFitness,
                                             final float maxFitness,
                                             final float avgFitness) {
    final float bestExpectedCopies = 2.0f;
    float a;
    float b;
    if (minFitness > (bestExpectedCopies*avgFitness - maxFitness) / (bestExpectedCopies - 1.0)) {
      final float fitnessRange = maxFitness - minFitness;
      a = (bestExpectedCopies - 1.0f) * avgFitness / fitnessRange;
      b = avgFitness * (maxFitness - bestExpectedCopies * avgFitness) / fitnessRange;
    } else {
      final float fitnessRange = avgFitness - minFitness;
      a = avgFitness / fitnessRange;
      b = -minFitness * avgFitness / fitnessRange;
    }
    
    return new float[] {a, b};
  }
  
  private static float scaleFitness(final float[] scalingCoefficients, final float fitness) {
    return scalingCoefficients[0] * fitness + scalingCoefficients[1];
  }
  
  private static float[] scaleFitnesses(final float[] fitnesses) {
    float sumOfFitnesses = fitnesses[0];
    float minFitness = fitnesses[0];
    float maxFitness = fitnesses[0];
    for (int index = 1; index < fitnesses.length; ++index) {
      if (fitnesses[index] < minFitness) {
        minFitness = fitnesses[index];
      } else if (fitnesses[index] > maxFitness) {
        maxFitness = fitnesses[index];
      }
      sumOfFitnesses += fitnesses[index];
    }
    final float avgFitness = sumOfFitnesses / fitnesses.length;

    final float[] scalingCoefficients = scalingCoefficients(minFitness, maxFitness, avgFitness);
    float[] scaledFitnesses = new float[fitnesses.length];
    for (int index = 0; index < fitnesses.length; ++index) {
      scaledFitnesses[index] = scaleFitness(scalingCoefficients, fitnesses[index]);
      System.out.println("Scaled Fitness " + index + ": " + scaledFitnesses[index]);
    }
    
    return scaledFitnesses;
  }

  /****************** STATIC HELPER METHODS ******************/

  public static float[][] loadExperimentalData(File file) {
    ArrayList<float[]> rowData = new ArrayList<float[]>();
    BufferedReader fin = null;
    try {
      fin = new BufferedReader(new FileReader(file));
      String seperator = "";
      String line = "";
      while ((line = fin.readLine()) != null) {
        // System.out.println(line);

        // determine the field seperator
        if (seperator.isEmpty()) {
          // System.out.print("Determining field seperator...");
          if (line.contains(",")) {
            seperator = ",";
          } else {
            seperator = " ";
          }
          // System.out.println(seperator);
        }
        
        // skip any header lines
        if (!line.isEmpty() && !line.matches("^[^0-9].*")) {
          // System.out.println("Not a header line.");
          String[] fields = line.split(seperator);
          int fieldIndex = 0;
          // skip the date field if it is present
          /*
          if (fields[fieldIndex].contains("/")) {
            System.out.println("Found date field.");
            ++fieldIndex;
          }*/
          // reformat the time field if necessary
          if (fields[fieldIndex].contains(".")) {
            /*
            System.out.println("Time field: " + fields[fieldIndex]);
            System.out.println("Time field . index: " + fields[fieldIndex].indexOf("\\."));
            System.out.println("Time field split size: " + fields[fieldIndex].split("\\.").length);
            */
            fields[fieldIndex] = fields[fieldIndex].split("\\.")[1];
          }
          float[] fieldData = new float[fields.length];
          // System.out.print("[");
          for (;fieldIndex<fields.length; ++fieldIndex) {
            fieldData[fieldIndex] = Float.parseFloat(fields[fieldIndex]);
            // System.out.print(fieldData[fieldIndex]);
            if (fieldIndex < (fields.length-1)) {
              // System.out.print(",");
            }
          }
          // System.out.println("]");
          rowData.add(fieldData);
        }
      }
    } catch (FileNotFoundException e) {
      System.out.println("ERROR: " + e.getMessage());
    } catch (IOException e) {
      System.out.println("ERROR: " + e.getMessage());
    } finally {
      if (fin != null) {
        try {
          fin.close();
        } catch (IOException e) {
          System.out.println("ERROR: " + e.getMessage());
        }
      }
    }
    
    // Convert the CSV rows into an array of signal data
    float[][] data = null;
    if (!rowData.isEmpty()) {
      data = new float[rowData.get(0).length][rowData.size()];
      Iterator<float[]> rowDatum = rowData.iterator();
      int sampleIndex = 0;
      float t0 = 0;
      while (rowDatum.hasNext()) {
        float[] fields = rowDatum.next();

        for (int signalIndex = 0; signalIndex<fields.length; ++signalIndex) {
          // convert us absolute times to s relative times
          if (signalIndex == 0) {
            if (sampleIndex == 0) {
              t0 = fields[signalIndex];
            }
            fields[signalIndex] = (fields[signalIndex] - t0) * (float) 1.0e-6;
          }
          data[signalIndex][sampleIndex] = fields[signalIndex];
          // System.out.println(data[signalIndex][sampleIndex]);
        }
        ++sampleIndex;
      }
      for (int signalIndex = 0; signalIndex<data.length; ++signalIndex) {
        data[signalIndex] = normalizeSignal(data[signalIndex]);
      }
    }

    return data;
  }
  
  /**
   *  Perform a normalized (subtract the mean and divide by the standard deviation)
   *  cross-correlation with zero lag.
   */
  public static float normalizedCrossCorrelation(float[] a, float[] b) {
    int length = a.length;
    if (b.length < a.length) {
      length = b.length;
    }
    // Calculate means
    float a_mean = 0;
    float b_mean = 0;
    for (int index=0; index<length; ++index) {
      a_mean += a[index];
      b_mean += b[index];
    }
    a_mean /= length;
    b_mean /= length;

    // Calculate standard deviations
    float[] a_diff = new float[length];
    float a_stddev = 0;
    float[] b_diff = new float[length];
    float b_stddev = 0;
    for (int index=0; index<length; ++index) {
      a_diff[index] = a[index] - a_mean;
      a_stddev += a_diff[index]*a_diff[index];

      b_diff[index] = b[index] - b_mean;
      b_stddev += b_diff[index]*b_diff[index];
    }
    a_stddev = (float) Math.sqrt(a_stddev);
    b_stddev = (float) Math.sqrt(b_stddev);

    // Calculate cross-correlation of normalized data only for zero lag
    float fitness = 0;
    for (int index=0; index<length; ++index) {
      fitness += a_diff[index] * b_diff[index];
    }
    if ((a_stddev == 0.0) || (b_stddev == 0.0)) {
      fitness = 0.0f;
    } else {
      fitness /= (a_stddev*b_stddev);
    }

    return fitness;
  }
  
  public static float rootMeanSquareError(float[] a, float[] b, int offset) {
    float ssr = 0.0f;
    for (int index = 0; index < a.length; ++index) {
      final float residual = a[index] - b[index+offset];
      ssr += residual*residual;
    }

    return Math.sqrt(ssr/a.length);
  }
  
  private static float[] signalStats(float[] signal) {
    float max = signal[0];
    float sum = signal[0];
    for (int index = 1; index < signal.length; ++index) {
      if (signal[index] > max) {
        max = signal[index];
      }
      sum += signal[index];
    }
    final float mean = sum / signal.length;

    return new float[] {max, mean};
  }
  
  private static float[] normalizeSignal(float[] signal) {
    final float[] signalStats = signalStats(signal);
    final float max = signalStats[0];
    final float mean = signalStats[1];
    float[] normalizedSignal = new float[signal.length];
    for (int index = 0; index < signal.length; ++index) {
      normalizedSignal[index] = (signal[index] - mean) / (max - mean);
    }
    
    return normalizedSignal;
  }
  
  public float[] evaluate(float[][] trueData) {
    // ### run the simulation batch job ###
    run();

    // ### play nice all of the simulations failed to produce data### 
    if (this.data == null) {
      float[] fitnesses = new float[this.parameterSets.size()];
      for (int index = 0; index < fitnesses.length; ++index) {
        fitnesses[index] = 1.0e-21f;
      }
      return fitnesses;
    }
    System.out.println("Length 1: " + this.data.length);
    System.out.println("Length 2: " + this.data[0].length);
    System.out.println("Length 3: " + this.data[0][0].length);

    // ### compare each set of plot data to the experimental data ###
    float sumOfFitnesses = 0.0f;
    float[] totalFitnesses = new float[this.data.length];
    for (int individualIndex = 0; individualIndex < this.data.length; ++individualIndex) {
      System.out.println("### Testing Individual " + individualIndex + " ###");
      final String[] genes = getGenes(this.parameterSets.get(individualIndex));
      for (int geneIndex = 0; geneIndex < genes.length; ++geneIndex) {
        System.out.print(genes[geneIndex] + " ");
      }
      System.out.println();
      totalFitnesses[individualIndex] = evaluateSimulationData(trueData, this.data[individualIndex],
                                                               this.t0s[individualIndex]);
      sumOfFitnesses += totalFitnesses[individualIndex];
    }

    // ### calculate the mean and standard deviation of the fitnesses for this generation ###
    final float mean = sumOfFitnesses / totalFitnesses.length;
    System.out.println("Mean Fitness: " + mean);
    float stddev = 0.0f;
    for (int individualIndex = 0; individualIndex < this.data.length; ++individualIndex) {
      stddev += Math.pow(totalFitnesses[individualIndex] - mean, 2);
    }
    stddev = (float) Math.sqrt(stddev / totalFitnesses.length);
    System.out.println("Fitness Standard Deviation: " + stddev);

    // ### convert fitnesses to values we can maximize ###
    for (int individualIndex = 0; individualIndex < this.data.length; ++individualIndex) {
      totalFitnesses[individualIndex] = (mean + 3 * stddev) - totalFitnesses[individualIndex];
      if (totalFitnesses[individualIndex] < 0) {
        totalFitnesses[individualIndex] = 0.0f;
      }
      System.out.println("Translated Fitness " + individualIndex + ": " + totalFitnesses[individualIndex]);
    }

    return totalFitnesses;
  }

  public static float evaluateSimulationData(float[][]      targetData,
                                             float[][]      simData,
                                             int            t0) {
    float totalFitness = 0.0f;
    System.out.println("Fitnesses: ");
    for (int plotIndex = 0; plotIndex < simData.length; ++plotIndex) {
    // for (int plotIndex = 0; plotIndex < 3; ++plotIndex) {
      // int plotIndex = 1;
      // Check if there was no plot and penalize if so
      float sampleSum = 0.0f;
      for (int index = 0; index < simData[plotIndex].length; ++index) {
        sampleSum += simData[plotIndex][index];
      }
      if (sampleSum == 0.0f) {
        totalFitness = 1.0e6f;
        break;
      }

      final float rootMeanSquareError = (float) rootMeanSquareError(
          simData[plotIndex],
          targetData[plotIndex],
          t0 / 2);  // 2us per sample
      System.out.print(rootMeanSquareError + " ");
      // System.out.print(individualFitness + " ");
      totalFitness += rootMeanSquareError;
    }
    System.out.println();
    System.out.println("Total Fitness: " + totalFitness);
    
    return totalFitness;
  }
    
  public static String[] getGenes(String paramsString) {
    return new String[] {
      paramsString.substring(0, 6),
      paramsString.substring(6,12),
      paramsString.substring(12,18),
      paramsString.substring(18,30),
      paramsString.substring(30,41),
      paramsString.substring(41,45),
      paramsString.substring(45,57),
      paramsString.substring(57,68),
      paramsString.substring(68,73),
      paramsString.substring(73,75)};
  }
  
  private static String repairDec(String dna) {
    final String[] genes = getGenes(dna);
    final String[][] minMaxGenes = new String[][]{
      getGenes(MIN_MAX_PARAM_STRINGS[0]), getGenes(MIN_MAX_PARAM_STRINGS[1])};
    // TODO: if gene is smaller, set to min. if gene is larger, set to max
    StringBuffer repairedDna = new StringBuffer(genes.length);
    boolean repaired = false;
    for (int geneIndex = 0; geneIndex < genes.length; ++geneIndex) {
      if (genes[geneIndex].compareTo(minMaxGenes[0][geneIndex]) < 0) {
        repairedDna.append(minMaxGenes[0][geneIndex]);
        repaired = true;
      } else if (genes[geneIndex].compareTo(minMaxGenes[1][geneIndex]) > 0) {
        repairedDna.append(minMaxGenes[1][geneIndex]);
        repaired = true;
      } else {
        repairedDna.append(genes[geneIndex]);
      }
    }
    if (repaired) {
      System.out.println("DNA Before: " + dna);
      System.out.println("Repaired DNA: " + repairedDna);
    }
    return repairedDna.toString();
  }

  private static String repair(String dna) {
    final int decDigitCount = dna.length() / 4;
    final String MAX_BINARY = "1001";  // decimal 9
    StringBuffer repairedDna = new StringBuffer(dna.length());
    for (int decDigitIndex = 0; decDigitIndex < decDigitCount; ++decDigitIndex) {
      final String decDigit = dna.substring(4*decDigitIndex, 4*decDigitIndex+4);
      String binQuad = null;
      if (decDigit.compareTo(MAX_BINARY) > 0) {
        binQuad = MAX_BINARY;
      } else {
        binQuad = decDigit;
      }
      repairedDna.append(binQuad);
    }
    return repairedDna.toString();
  }

  private static String generateRandomParameter(final String templateGene,
                                                final long minValue,
                                                final long maxValue) {
    Random rand = new Random();
    final long range = maxValue - minValue;
    // System.out.println("Range: " + range);
    final long parameter = (long) Math.round(range *rand.nextDouble()) + minValue;
    // System.out.println("Parameter " + paramIndex + ": " + parameter);
    final String paramString = String.valueOf(parameter);
    /*
    System.out.println("Param String: " + paramString);
    System.out.println("Gene String: " + genes[paramIndex]);
    */
    final int lengthDelta = templateGene.length() - paramString.length();
    final String padding = templateGene.substring(0, lengthDelta);

    return padding + paramString;
  }

  public static String generateRandomIndividual() {
    ArrayList<Long> parameters = new ArrayList<Long>(PARAMETER_COUNT);
    StringBuffer dna = new StringBuffer();
    for (int paramIndex = 0; paramIndex < PARAMETER_COUNT; ++paramIndex) {
      dna.append(generateRandomParameter(PARAM_STRING_TEMPLATES[paramIndex],
                                         MIN_MAX_PARAMS[0][paramIndex],
                                         MIN_MAX_PARAMS[1][paramIndex]));
    }
    
    return dna.toString();
  }

  public static String mutateDec(String dna, double chancePerBase) {
    Random rand = new Random();
    StringBuffer mutatedDna = new StringBuffer(dna.length());
    // for (int index=0; index<dna.length(); ++index) {
    for (int charIndex=0; charIndex<dna.length(); ++charIndex) {
      if (rand.nextFloat() <= chancePerBase) {
        if ((charIndex == 0) || (charIndex == 6)) {
          // 10s place in r_S0 and r_S5
          // get random digit between 0 and 1
          mutatedDna.append((int) Math.round(rand.nextDouble()));
        } else if (charIndex == 12) {
          // 10s place in z_S1
          // digit must be 0
          mutatedDna.append(0);
        } else if (charIndex == 13) {
          //100s place in z_S1
          // get random digit between 0 and 8
          mutatedDna.append((int) Math.round(8 * rand.nextDouble()));
          /*
          char digit = dna.charAt(charIndex);
          if (((int) Math.round(rand.nextDouble())) > 0) {
            if (digit == '8') {
              mutatedDna.append(0);
            } else {
              mutatedDna.append(digit+1);
            }
          } else {
            if (digit == '0') {
              mutatedDna.append(8);
            } else {
              mutatedDna.append(digit-1);
            }
          }
          */
        } else {
          // get random digit between 0 and 9
          mutatedDna.append((int) Math.round(9 * rand.nextDouble()));
          // increment/decrement
          /*
          char digit = dna.charAt(charIndex);
          if (((int) Math.round(rand.nextDouble())) > 0) {
            if (digit == '9') {
              mutatedDna.append(0);
            } else {
              mutatedDna.append(digit+1);
            }
          } else {
            if (digit == '0') {
              mutatedDna.append(9);
            } else {
              mutatedDna.append(digit-1);
            }
          }
          */
        }
      } else {
        mutatedDna.append(dna.charAt(charIndex));
      }
    }
    
    return repairDec(mutatedDna.toString());
  }

  public static String mutate(String dna, double chancePerBase) {
    Random rand = new Random();
    String binaryDna = binaryEncodeParamsString(dna);
    StringBuffer mutatedDna = new StringBuffer(dna.length());
    for (int bitIndex=0; bitIndex<binaryDna.length(); ++bitIndex) {
      char currentBit = binaryDna.charAt(bitIndex);
      if (rand.nextFloat() <= chancePerBase) {
        if (currentBit == '0') {
          mutatedDna.append((char) '1');
        } else {
          mutatedDna.append((char) '0');
        }
      } else {
        mutatedDna.append((char) currentBit);
      }
    }
    
    return repairDec(binaryDecodeParamsString(repair(mutatedDna.toString())));
  }
  
  public static String[] crossover(String dna1, String dna2) {
    // System.out.println("Parent DNA: " + dna1 + " " + dna2);
    String[] encodedDna = new String[] {
      binaryEncodeParamsString(dna1), binaryEncodeParamsString(dna2)};
    // System.out.println("dna1: " + encodedDna[0] + "\n" + "dna2: " + encodedDna[1]);
    Random rand = new Random();
    String[] crossedDna = null;
    // 60% probability of crossing
    if (rand.nextFloat() < 0.6f) {
      int dnaLength = encodedDna[0].length();
      int crossSiteIndex = (int) Math.round((dnaLength-2) * rand.nextFloat()) + 1;
      // System.out.println("Cross Site Index: " + crossSiteIndex);
      crossedDna = new String[] {
        repairDec(binaryDecodeParamsString(repair(
          encodedDna[0].substring(0, crossSiteIndex) + encodedDna[1].substring(crossSiteIndex)))),
        repairDec(binaryDecodeParamsString(repair(
          encodedDna[1].substring(0, crossSiteIndex) + encodedDna[0].substring(crossSiteIndex))))};
    } else {
      crossedDna = new String[] {dna1, dna2};
    }

    // System.out.println("Crossed DNA: " + crossedDna[0] + " " + crossedDna[1]);
    return crossedDna;
  }
  
  public static String mate(String dna1, String dna2) {
    /*
    int[] startIndicies = new int[] {0, 6, 12, 18, 30, 41, 45, 56, 67};
    int[] endIndicies = new int[] {6, 12, 18, 30, 41, 45, 56, 67, 71};
    */
    String[] parentDna = new String[] {dna1, dna2};
    Random rand = new Random();
    StringBuffer dna = new StringBuffer(dna1.length());
    for (int baseIndex = 0; baseIndex < dna1.length(); ++baseIndex) {
      int parentIndex = (int) Math.round(rand.nextDouble());
      // dna.append(parentDna[parentIndex].substring(startIndicies[geneIndex], endIndicies[geneIndex]));
      dna.append(parentDna[parentIndex].substring(baseIndex, baseIndex+1));
    }

    return repair(dna.toString());
  }

  // From https://www.ee.columbia.edu/~ronw/code/dev/MEAPsoft/src/com/meapsoft/DSP.java
  /**
   * Computes the cross correlation between sequences a and b.
   * maxlag is the maximum lag to
   */
   /*
  public static double[] xcorr(double[] a, double[] b, int maxlag)
  {
    double[] y = new double[2*maxlag+1];
    Arrays.fill(y, 0);
    
    for(int lag = b.length-1, idx = maxlag-b.length+1; 
      lag > -a.length; lag--, idx++)
    {
      if(idx < 0)
        continue;
      
      if(idx >= y.length)
        break;

      // where do the two signals overlap?
      int start = 0;
      // we can't start past the left end of b
      if(lag < 0) 
      {
        //System.out.println("b");
        start = -lag;
      }

      int end = a.length-1;
      // we can't go past the right end of b
      if(end > b.length-lag-1)
      {
        end = b.length-lag-1;
        //System.out.println("a "+end);
      }

      //System.out.println("lag = " + lag +": "+ start+" to " + end+"   idx = "+idx);
      for(int n = start; n <= end; n++)
      {
        //System.out.println("  bi = " + (lag+n) + ", ai = " + n); 
        y[idx] += a[n]*b[lag+n];
      }
      //System.out.println(y[idx]);
    }

    return(y);
  }
  */

  /*
  public double evaluate(org.jgap.IChromosome solution) {
    // Extract alleles from solution
    // Set simulation parameters from alleles
    // Run the simulation
    // Return fitness() + 1
    return 0.0;
  }
*/

  /****************** Test Functions ******************/

  public static String[] getTestGenes(String paramsString) {
    return new String[] {
      paramsString.substring(0, 36),
      paramsString.substring(36)};
  }

  public float[] getTestParams(int paramSetIndex) {
    String paramsString = this.parameterSets.get(paramSetIndex);
    // System.out.println("Params String: " + paramsString);
    String[] genes = getTestGenes(paramsString);
    // System.out.println("Genes: " + genes[0] + " " + genes[1]);
    float[] parameters = new float[genes.length];
    for (int paramIndex = 0; paramIndex < parameters.length; ++paramIndex) {
      parameters[paramIndex] = Float.parseFloat(genes[paramIndex].substring(0, 1) + "." +
                                                genes[paramIndex].substring(1))
                             -5.0f;
    }

    return parameters;
  }

  public float[] testEvaluate() {
    testRun();

    // final float MIN_FITNESS = -12.64241f;  // for 2D Ackley's function
    final float MIN_FITNESS = -40016.0f;  // for the 2D Rosenbrock function
    // final float MIN_FITNESS = -45.0f;  // for f(x,y) = 5-x^2-y^2
    float[] fitnesses = new float[this.data.length];
    for (int individualIndex = 0; individualIndex < this.data.length; ++individualIndex) {
      float x = this.data[individualIndex][0][0];
      float y = this.data[individualIndex][0][1];
      System.out.println(individualIndex + "(x, y) = (" + x + "," + y + ")");

      // calculate the value of 2D Ackley's function -- max at (x, y) = (0, 0)
      /*
      fitnesses[individualIndex] = -(-20.0f * (float) Math.exp(-0.2*Math.sqrt(0.5*(x*x+y*y)))
                                 - (float) Math.exp(0.5*(Math.cos(2.0*Math.PI*x) + Math.cos(2.0*Math.PI*y)))
                                 + (float) Math.exp(1)
                                 + 20.0f) - MIN_FITNESS;
                                 */
      // calculate the value of the 2D Rosenbrock function -- max at (x, y) = (1, 1)
      fitnesses[individualIndex] = -(100.0f * (y - x*x)*(y - x*x) + (1.0f-x)*(1.0f-x)) - MIN_FITNESS;
      // calcuate the value of a simple, inverted, 2D parabola -- max at (x, y) = (0, 0)
      // fitnesses[individualIndex] = (float) (5-x*x-y*y - MIN_FITNESS);
    }
    float[] scaledFitnesses = scaleFitnesses(fitnesses);
  
    return scaledFitnesses;
  }
  
  public static String testMate(String dna1, String dna2) {
    String[] parentDna = new String[] {dna1, dna2};
    Random rand = new Random();
    StringBuffer dna = new StringBuffer(dna1.length());
    for (int baseIndex = 0; baseIndex < dna1.length(); ++baseIndex) {
      int parentIndex = (int) Math.round(rand.nextDouble());
      dna.append(parentDna[parentIndex].substring(baseIndex, baseIndex+1));
    }

    return dna.toString();
  }

  public static String testMutate(String dna, double chancePerBase) {
    Random rand = new Random();
    StringBuffer mutatedDna = new StringBuffer(dna.length());
    for (int charIndex=0; charIndex<dna.length(); ++charIndex) {
      if (rand.nextFloat() <= chancePerBase) {
        // get random digit between 0 and 9
        mutatedDna.append((int) Math.round(9 * rand.nextDouble()));
      } else {
        mutatedDna.append(dna.charAt(charIndex));
      }
    }
    
    return mutatedDna.toString();
  }

  public void testRun() {
    final int SOLUTION_COUNT = parameterSets.size();
    this.data = new float[SOLUTION_COUNT][1][2];
    for (int solutionIndex = 0; solutionIndex < SOLUTION_COUNT; ++solutionIndex) {
      // convert the parameter string to x and y floats
      float[] parameters = getTestParams(solutionIndex);
      for (int paramIndex = 0; paramIndex < 2; ++paramIndex) {
        this.data[solutionIndex][0][paramIndex] = parameters[paramIndex];
      }
    }
  }

  public static void testMain(String[] args) {
    final long DURATION = 200;  // us
    final int POPULATION_SIZE = 500;
    final int GENERATION_COUNT = 9000000;  // Good value: 100*# variables = 900
    final float MUTATION_RATE = 1.0f/POPULATION_SIZE/2.0f;  // arbitrarily taken from Goldberg, Ch. 1
    // final float MUTATION_RATE = 0.033f;

    System.out.println(binaryEncodeParamsString(MIN_MAX_PARAM_STRINGS[0]));
    System.out.println(binaryEncodeParamsString(MIN_MAX_PARAM_STRINGS[1]));
    // #### Load or create the model ####
    HprfRfHammer2DGeneticSolver sim = new HprfRfHammer2DGeneticSolver();

    // ### Create a completely random initial population ###
    sim.clearParams();
    final String template = "00000000000000000000000000000000000000000000000000000000000000000000000";
    for (int individualIndex = 0; individualIndex < POPULATION_SIZE; ++individualIndex) {
      // final String paramString = HprfRfHammer2DGeneticSolver.mutate(template, 0.5);
      final String paramString = HprfRfHammer2DGeneticSolver.mutateDec(template, 1.0);
      System.out.println("Base-10 Param String:\n" + paramString);
      final String encodedParamString = binaryEncodeParamsString(paramString);
      System.out.println("Base-2 Param String:\n" + encodedParamString);
      final String decodedParamString = binaryDecodeParamsString(encodedParamString);
      System.out.println("Decoded Param String:\n" + decodedParamString);
      sim.addParamsFromString(paramString);
    }

    float[] fitnesses = sim.testEvaluate();
    /*
    System.out.println("### Initial Population Fitnesses ###");
    for (int fitIndex = 0; fitIndex < fitnesses.length; ++fitIndex) {
      System.out.println("Individual " + fitIndex + ": " + fitnesses[fitIndex]);
      System.out.println(sim.getParamsAsString(fitIndex));
    }
    */

    for (int generationIndex = 0; generationIndex < GENERATION_COUNT; ++generationIndex) {
      // System.out.println("# Fitnesses: " + fitnesses.length);
      // sum the fitnesses
      float sumOfFitnesses = 0.0f;
      for (int fitIndex = 0; fitIndex < fitnesses.length; ++fitIndex) {
        sumOfFitnesses += fitnesses[fitIndex];
      }
      final float meanFitness = sumOfFitnesses / fitnesses.length;
      System.out.println("Mean Scaled Fitness: " + meanFitness);
      float stddevOfFitnesses = 0.0f;
      for (int fitIndex = 0; fitIndex < fitnesses.length; ++fitIndex) {
        stddevOfFitnesses += Math.pow((fitnesses[fitIndex] - meanFitness), 2);
      }
      stddevOfFitnesses = (float) Math.sqrt(stddevOfFitnesses / fitnesses.length);
      System.out.println("Std. Dev.: " + stddevOfFitnesses);
      
      // select the fittest individuals from the current generation
      ArrayList<String> selections = sim.selectFittest(fitnesses);
      System.out.println("Selection Count: " + selections.size());
      // pair off and mate the fittest individuals to create the next generation
      ArrayList<String> population = new ArrayList<String>(fitnesses.length);
      final Random rand = new Random();
      while (selections.size() > 0) {
        String[] children = HprfRfHammer2DGeneticSolver.crossover(
          selections.remove((int) Math.round((selections.size() - 1) * rand.nextFloat())),
          selections.remove((int) Math.round((selections.size() - 1) * rand.nextFloat())));
        for (int childIndex = 0; childIndex < 2; ++childIndex) {
          population.add(HprfRfHammer2DGeneticSolver.mutate(children[childIndex], MUTATION_RATE));
        }
      }

      // add the next generation to the simulation
      sim.clearParams();
      Iterator<String> iter = population.iterator();
      while (iter.hasNext()) {
        sim.addParamsFromString(iter.next());
      }

      // run the simulation and generate a fitness for each individual
      fitnesses = sim.testEvaluate();

      System.out.println("### Generation " + generationIndex + " Fitnesses ###");
      int maxFitnessIndex = 0;
      for (int fitIndex = 0; fitIndex < fitnesses.length; ++fitIndex) {
        if (fitnesses[fitIndex] > fitnesses[maxFitnessIndex]) {
          maxFitnessIndex = fitIndex;
        }
        /*
        System.out.println("Individual " + fitIndex + ": " + fitnesses[fitIndex]);
        System.out.println(sim.getParamsAsString(fitIndex));
        */
      }
      System.out.println("Fittest Individual: " + maxFitnessIndex + " (" + fitnesses[maxFitnessIndex] + ")");
      final String[] genes = getGenes(sim.getParamsAsString(maxFitnessIndex));
      for (int geneIndex = 0; geneIndex < genes.length; ++geneIndex) {
        System.out.print(genes[geneIndex] + " ");
      }
      System.out.println();
    }  // generation loop
  }

  public static float[][] getTestData() {
    /* Generated on July 17, 2015
     * Duration: 500us
     * Param String: 10000010000000000014296543266677666555444798711011011014247500000000889212
     */
    return new float[][] {
     {-0.008742699f,-0.008742699f,-0.008740917f,-0.008732894f,-0.008716586f,-0.008724117f,-0.008868744f,-0.009297008f,-0.009869774f,-0.009617784f,-0.006570394f,0.0011469272f,0.012857251f,0.023203356f,0.023088302f,0.005038482f,-0.029080689f,-0.06291971f,-0.067871526f,-0.014080698f,0.11345681f,0.3017586f,0.5084572f,0.67595994f,0.7558853f,0.730926f,0.6212814f,0.47129014f,0.32358605f,0.19558643f,0.07202603f,-0.08172123f,-0.29283983f,-0.5565789f,-0.82941836f,-1.0453224f,-1.1466974f,-1.112581f,-0.96794844f,-0.76888645f,-0.57216454f,-0.40659335f,-0.26256096f,-0.105601534f,0.094478644f,0.33696336f,0.58260393f,0.7695976f,0.8432111f,0.7815136f,0.6034191f,0.35672823f,0.09508791f,-0.1427221f,-0.34125835f,-0.50280935f,-0.6336157f,-0.731401f,-0.7818862f,-0.76607007f,-0.67338884f,-0.5121278f,-0.3099086f,-0.10298985f,0.07996273f,0.2302711f,0.36030003f,0.49209905f,0.6397779f,0.7955312f,0.92738014f,0.9903784f,0.9458762f,0.77936393f,0.5084228f,0.17748864f,-0.1574871f,-0.44710284f,-0.6615432f,-0.7940573f,-0.8551868f,-0.86205006f,-0.8284223f,-0.7599607f,-0.65571564f,-0.5138281f,-0.33755407f,-0.1381806f,0.06647268f,0.25671852f,0.4160141f,0.5336964f,0.60451806f,0.62615556f,0.59704757f,0.5165064f,0.38733506f,0.21930034f,0.030951537f,-0.15201229f,-0.30179247f,-0.39507362f,-0.41846973f,-0.37081888f,-0.26193208f,-0.108986326f,0.06754356f,0.24693966f,0.4099011f,0.53942966f,0.6218705f,0.6485641f,0.61770934f,0.5354005f,0.41481733f,0.27325127f,0.12768874f,-0.009501605f,-0.1330614f,-0.24442442f,-0.3484795f,-0.4488361f,-0.5436345f,-0.6238037f,-0.6746917f,-0.6805758f,-0.63032144f,-0.521942f,-0.3642491f,-0.17498048f,0.023796413f,0.21118128f,0.3718083f,0.49768683f,0.5872638f,0.642579f,0.66602045f,0.6581254f,0.61722875f,0.54088104f,0.42823768f,0.28233021f,0.11131408f,-0.071718484f,-0.25022304f,-0.40647712f,-0.5240392f,-0.5899193f,-0.5962102f,-0.5411032f,-0.4292738f,-0.27159762f,-0.08413781f,0.1135977f,0.30097646f,0.45879897f,0.5716294f,0.6294555f,0.6284081f,0.57052326f,0.46276838f,0.3156684f,0.14185466f,-0.04525916f,-0.23259045f,-0.4082804f,-0.5623547f,-0.6871109f,-0.777191f,-0.82935953f,-0.842118f,-0.815371f,-0.7503558f,-0.6499214f,-0.5190362f,-0.36521867f,-0.19852436f,-0.030839022f,0.12550882f,0.25950268f,0.36305195f,0.43195936f,0.4659538f,0.46783715f,0.44211784f,0.39366832f,0.3268929f,0.24564683f,0.15382011f,0.05622207f,-0.040715374f,-0.12882476f,-0.19886325f,-0.24184722f,-0.25072873f,-0.22197026f,-0.15660816f,-0.0605465f,0.05597367f,0.17959027f,0.29545698f,0.38899508f,0.4474936f,0.46138373f,0.42511296f,0.33763024f,0.20252857f,0.027890386f,-0.17415203f,-0.38816056f,-0.5964674f,-0.780579f,-0.9228381f,-1.0081933f,-1.0258616f,-0.97065115f,-0.8437392f,-0.65278804f,-0.4113794f,-0.13784677f,0.14634643f,0.41850922f,0.6567325f,0.8418478f,0.95918757f,1.0f,0.9623656f,0.85147476f,0.6791671f,0.46271828f,0.2229819f,-0.01788447f,-0.23878047f,-0.42175147f,-0.5534967f,-0.6261161f,-0.637089f,-0.58864766f,-0.48682347f,-0.3404703f,-0.16049728f,0.04060142f,0.24899651f,0.44976094f,0.6273436f,0.7664673f,0.85348743f,0.8780567f,0.8347779f,0.72444886f,0.5545389f,0.33869463f,0.095289044f,-0.15474147f,-0.39036778f,-0.5927884f,-0.74721044f},
     {-0.009115097f,-0.009115097f,-0.00944423f,-0.010972984f,-0.015239317f,-0.023002222f,-0.032117106f,-0.0376816f,-0.035581253f,-0.025555082f,-0.010063964f,0.008184664f,0.025059713f,0.032102928f,0.018892916f,-0.017706994f,-0.06453566f,-0.091815464f,-0.063785195f,0.044012085f,0.22815908f,0.4501324f,0.6462906f,0.75554645f,0.75142854f,0.6568694f,0.52854705f,0.41891375f,0.34108785f,0.26075548f,0.120939605f,-0.11671584f,-0.44122118f,-0.79140466f,-1.0840149f,-1.2517178f,-1.2685394f,-1.1519499f,-0.94633996f,-0.7006199f,-0.45090353f,-0.21343462f,0.0118312575f,0.22993368f,0.437951f,0.61786854f,0.74001205f,0.77541065f,0.7097859f,0.55105764f,0.32654494f,0.072103225f,-0.18042834f,-0.41070658f,-0.6087585f,-0.76704156f,-0.87219816f,-0.9039274f,-0.84394646f,-0.6903415f,-0.46702242f,-0.2193773f,0.0042712633f,0.17557633f,0.29910806f,0.40608612f,0.53151315f,0.6881183f,0.8529796f,0.9751449f,1.0f,0.8965811f,0.6731345f,0.3737446f,0.05949732f,-0.21505356f,-0.41882795f,-0.54851174f,-0.6211219f,-0.6605811f,-0.6850221f,-0.6994923f,-0.6955835f,-0.65685135f,-0.56733024f,-0.41993108f,-0.22181836f,0.00507466f,0.2286021f,0.41447973f,0.53510034f,0.5761581f,0.5391145f,0.43915957f,0.2998796f,0.1467479f,0.0015822347f,-0.120576434f,-0.21176007f,-0.2691447f,-0.29215953f,-0.28005832f,-0.23103519f,-0.14326087f,-0.017391203f,0.14056112f,0.31736067f,0.49280944f,0.64255226f,0.7429167f,0.7764442f,0.7363603f,0.6284421f,0.46954414f,0.2831656f,0.093473524f,-0.080278456f,-0.22717185f,-0.3449451f,-0.43741968f,-0.51026684f,-0.5668845f,-0.6060024f,-0.62186617f,-0.60678446f,-0.55486506f,-0.4653122f,-0.34388262f,-0.20189944f,-0.053269193f,0.08923366f,0.21689202f,0.32607037f,0.4170937f,0.4916592f,0.5500341f,0.58929485f,0.60334766f,0.58469504f,0.5272314f,0.4290095f,0.29401365f,0.1323967f,-0.040814046f,-0.20807056f,-0.35201228f,-0.45785326f,-0.51501614f,-0.5179138f,-0.4660654f,-0.3638676f,-0.22027628f,-0.048445523f,0.13486557f,0.31042773f,0.4584363f,0.56109065f,0.6054091f,0.5855579f,0.5039627f,0.37076533f,0.20169717f,0.014954879f,-0.17202275f,-0.34500206f,-0.49408808f,-0.6135115f,-0.7004964f,-0.7538819f,-0.7731278f,-0.7580612f,-0.7093535f,-0.62940246f,-0.5231414f,-0.39833623f,-0.2651146f,-0.13473038f,-0.017813269f,0.07747543f,0.14689784f,0.1907892f,0.21366794f,0.22283822f,0.22630323f,0.23053494f,0.23873042f,0.2500823f,0.26030684f,0.2632894f,0.25334743f,0.2273993f,0.18634988f,0.13525222f,0.08219744f,0.036282893f,0.0052769077f,-0.0063448236f,0.0015262245f,0.024664916f,0.055544417f,0.08500571f,0.10409314f,0.10562318f,0.08519908f,0.041575678f,-0.023548925f,-0.106123656f,-0.20048538f,-0.30004805f,-0.39781523f,-0.48669133f,-0.5596624f,-0.6099767f,-0.63146716f,-0.6191032f,-0.56976116f,-0.48308322f,-0.36219865f,-0.21404997f,-0.049118698f,0.11952456f,0.2777454f,0.41218388f,0.5120437f,0.57047015f,0.58525115f,0.5587457f,0.49712393f,0.4091511f,0.3048129f,0.19406293f,0.08588242f,-0.012274405f,-0.09468675f,-0.15733051f,-0.19772007f,-0.2147255f,-0.2084165f,-0.17994255f,-0.13143454f,-0.06591154f,0.0128131425f,0.10022135f,0.19116989f,0.27996933f,0.36049312f,0.42636606f,0.47128472f,0.48950118f,0.4764583f,0.42950964f,0.34860465f,0.236791f,0.10038692f,-0.051279154f,-0.20658925f,-0.35275468f},
     {0.0029016118f,0.0029016118f,0.002900201f,0.0028929182f,0.0028731832f,0.0028547284f,0.0029245054f,0.0032752296f,0.004100143f,0.005310879f,0.0063115777f,0.0062700068f,0.005112891f,0.0048481603f,0.010226902f,0.027682228f,0.06215019f,0.112725146f,0.16958182f,0.21506472f,0.2302429f,0.20482835f,0.1452053f,0.07496665f,0.025798403f,0.0221618f,0.06761258f,0.14074144f,0.20392922f,0.22081733f,0.17307398f,0.06738811f,-0.07032405f,-0.20883512f,-0.32492033f,-0.40927616f,-0.46405274f,-0.4961091f,-0.51122344f,-0.51182806f,-0.49733663f,-0.46479097f,-0.4090342f,-0.32384583f,-0.2056472f,-0.058762047f,0.10178449f,0.25380418f,0.37737083f,0.46427047f,0.52005947f,0.555908f,0.5752261f,0.5650015f,0.49947315f,0.35543975f,0.13014139f,-0.14957333f,-0.4326231f,-0.6604009f,-0.7854202f,-0.7830265f,-0.65316516f,-0.41545072f,-0.10304502f,0.24151507f,0.5694436f,0.83073306f,0.982655f,1.0f,0.88286024f,0.6579196f,0.37177187f,0.078110896f,-0.17694558f,-0.36644724f,-0.48720723f,-0.5558418f,-0.5977485f,-0.63339126f,-0.6677422f,-0.6877752f,-0.66943157f,-0.5908724f,-0.4454811f,-0.24802232f,-0.03089146f,0.16738541f,0.31531933f,0.39687353f,0.4120064f,0.37097672f,0.28725475f,0.17368586f,0.0433777f,-0.08703935f,-0.19598177f,-0.26022002f,-0.26235026f,-0.19857676f,-0.082067445f,0.06073237f,0.19961536f,0.3115312f,0.38774484f,0.4341309f,0.46505824f,0.49435166f,0.52778757f,0.56025326f,0.57808524f,0.564722f,0.50679713f,0.39834133f,0.24222308f,0.049347572f,-0.16323683f,-0.37441176f,-0.56102043f,-0.7001995f,-0.7724001f,-0.7650104f,-0.6758118f,-0.5150544f,-0.30502772f,-0.076666094f,0.1362813f,0.30406877f,0.40686023f,0.43802953f,0.4039751f,0.3206964f,0.20840447f,0.0859963f,-0.03286932f,-0.14104488f,-0.23646803f,-0.3189443f,-0.3866235f,-0.43368062f,-0.45021012f,-0.42451045f,-0.34706447f,-0.21489772f,-0.03481118f,0.17571232f,0.39044756f,0.5788906f,0.71218365f,0.7690208f,0.7400701f,0.6298653f,0.45580173f,0.2445954f,0.027131412f,-0.16708842f,-0.3147299f,-0.4017759f,-0.4246322f,-0.38926408f,-0.30883896f,-0.20058982f,-0.08265878f,0.028463885f,0.11956797f,0.1814046f,0.20866892f,0.19965507f,0.15590309f,0.08197858f,-0.0146889f,-0.12412412f,-0.23449728f,-0.33326954f,-0.4086751f,-0.4512625f,-0.4551697f,-0.41889f,-0.34543937f,-0.24198176f,-0.11903741f,0.010617235f,0.13330506f,0.2360288f,0.30818665f,0.34320512f,0.33974236f,0.3020801f,0.2394429f,0.1642502f,0.08963518f,0.026828289f,-0.016914165f,-0.039296437f,-0.042670697f,-0.032604944f,-0.015837377f,0.0017272652f,0.016389469f,0.027049601f,0.034612533f,0.040700585f,0.04617341f,0.049993604f,0.048856813f,0.03776424f,0.011422007f,-0.033915233f,-0.09870403f,-0.17903344f,-0.2666081f,-0.34990135f,-0.4162437f,-0.45432258f,-0.4564342f,-0.4198986f,-0.34730634f,-0.24562615f,-0.12453672f,0.005456664f,0.13457021f,0.2546018f,0.35905176f,0.44273657f,0.5012882f,0.53091896f,0.5286651f,0.4930823f,0.42515114f,0.3290389f,0.21238594f,0.08592115f,-0.03760162f,-0.14497802f,-0.22435468f,-0.26693848f,-0.26832417f,-0.22922742f,-0.15550433f,-0.0574407f,0.051603153f,0.15700777f,0.24490091f,0.30403528f,0.3272228f,0.3120834f,0.26100165f,0.18035215f,0.079211414f,-0.03212775f,-0.14351268f,-0.24591348f,-0.33193386f,-0.39591426f,-0.43381542f,-0.4431128f,-0.42288157f}};  }

  public static void testFitness(String[] args) {
    float[][] testData = getTestData();
    float fitness = 0.0f;
    int t0 = 0;

    for (int i = 0; i < testData.length; ++i) {
      final float[] stats = signalStats(testData[i]);
      System.out.println("Test Signal " + i + " max value: " + stats[0]);
      System.out.println("Test Signal " + i + " mean value: " + stats[1]);
    }

    fitness = evaluateSimulationData(testData, testData, t0);
    System.out.println("Fitness for test data: " + fitness);

    float[][] negatedData = new float[testData.length][testData[0].length];
    for (int i = 0; i < testData.length; ++i) {
      for (int j = 0; j < testData.length; ++j) {
        negatedData[i][j] = -testData[i][j];
      }
    }
    fitness = evaluateSimulationData(testData, negatedData, t0);
    System.out.println("Fitness for negated test data: " + fitness);

    /*
    float[][] negatedData = new float[testData.length][testData[0].length];
    for (int i = 0; i < testData.length; ++i) {
      for (int j = 0; j < testData.length; ++j) {
        negatedData[i][j] = -testData[i][j];
      }
    }
    */
  }

  /****************** MAIN ******************/
  public static void main0(String[] args) {
    if ((args == null) || (args.length != 2)) {
      System.out.println("Usage: comsolrun HprfRfHammer2DGeneticSolver.java <POPULATION_SIZE> <GENERATION_COUNT>");
      return;
    }
    final int POPULATION_SIZE = Integer.parseInt(args[0]);
    final int GENERATION_COUNT = Integer.parseInt(args[1]);
    final int DURATION = 200;  // us
    // final float MUTATION_RATE = 1.0f/POPULATION_SIZE;  // arbitrarily taken from Goldberg, Ch. 1
    final float MUTATION_RATE = 0.001f;
    final float RELATIVE_TOLERANCE = 1.0e-6f;

    final boolean CREATE_NEW_MODEL = true;
    final boolean LOAD_TEST_DATA = false;

    // #### Load the real data ####
    float[][] trueData = null;
    if (LOAD_TEST_DATA) {
      trueData = getTestData();
    } else {
      File labviewDataFile
        = new File("C:\\Users\\plane\\Desktop\\COMSOL\\Analysis\\RF Hammer",
                   "RF_hammer_36MV-M_500us.csv");
      // S0, S1, and S5 are at indicies 1, 2, and 6 respectively (index 0 is time)
      System.out.print("Loading experimental data...");
      float[][] labviewData = HprfRfHammer2DGeneticSolver.loadExperimentalData(labviewDataFile);
      System.out.println("Done");
      // time is placed in the first row, so S0 has index 1 etc...
      final int[] labviewIndicies = new int[] {1, 6, 2};
      System.out.println("Experimental Data:");
      for (int signalIndex = 0; signalIndex < 3; ++signalIndex) {
        System.out.print("signal_" + signalIndex + " = np.array((");
        for (int sampleIndex = 0; sampleIndex < labviewData[signalIndex].length; ++sampleIndex) {
          if (sampleIndex > 0) {
            System.out.print(",");
          }
          System.out.print(labviewData[labviewIndicies[signalIndex]][sampleIndex]);
        }
        System.out.println("))");
      }
      trueData = new float[][] {
        labviewData[labviewIndicies[0]],
        labviewData[labviewIndicies[1]],
        labviewData[labviewIndicies[2]]};
    }
        
    // #### Load or create the model ####
    File modelFile = new File(
      "C:\\msys64\\home\\plane\\Development\\COMSOL\\HprfRfHammer2DGeneticSolver_Model0.mph");
    HprfRfHammer2DGeneticSolver sim = null;
    if (CREATE_NEW_MODEL) {
      sim = HprfRfHammer2DGeneticSolver.create();
      sim.setDuration(DURATION);
      sim.build();
    } else {
      try {
        // File modelFile = sim.save();

        sim = HprfRfHammer2DGeneticSolver.load(modelFile);
      } catch(IOException ioe) {
        System.out.println("ERROR: " + ioe.getMessage());
        System.out.flush();
        ioe.printStackTrace();
        try {
            Thread.sleep(1000);                 //1000 milliseconds is one second.
        } catch(InterruptedException ex) {
            Thread.currentThread().interrupt();
        }
      }
      sim.setDuration(DURATION);
    }

    /* Training Parameters
      "10000010000000000013400000000075000000000807099000000000440000000008990"
    */

    // ### Create a completely random initial population ###
    sim.clearParams();
    for (int individualIndex = 0; individualIndex < POPULATION_SIZE; ++individualIndex) {
      final String paramString = generateRandomIndividual();
      // final String paramString = HprfRfHammer2DGeneticSolver.mutate(template, 0.5);
      // final String paramString = HprfRfHammer2DGeneticSolver.mutateDec(template, 1.0);
      System.out.println("Base-10 Param String:\n" + paramString);
      final String encodedParamString = binaryEncodeParamsString(paramString);
      System.out.println("Base-2 Param String:\n" + encodedParamString);
      final String decodedParamString = binaryDecodeParamsString(encodedParamString);
      System.out.println("Decoded Param String:\n" + decodedParamString);
      sim.addParamsFromString(paramString);
    }

    float[] fitnesses = sim.evaluate(trueData);
    float[] scaledFitnesses = scaleFitnesses(fitnesses);
    float avgFitness = 0.0f;
    System.out.println("### Initial Population Fitnesses ###");
    for (int fitIndex = 0; fitIndex < fitnesses.length; ++fitIndex) {
      avgFitness += fitnesses[fitIndex];
      System.out.println("Individual " + fitIndex + ": " + scaledFitnesses[fitIndex]);
      System.out.println(sim.getParamsAsString(fitIndex));
    }
    avgFitness /= fitnesses.length;

    int maxFitnessIndex = 0;

    for (int generationIndex = 0; generationIndex < GENERATION_COUNT; ++generationIndex) {
      // System.out.println("# Fitnesses: " + scaledFitnesses.length);
      // sum the fitnesses
      float sumOfFitnesses = 0.0f;
      for (int fitIndex = 0; fitIndex < scaledFitnesses.length; ++fitIndex) {
        sumOfFitnesses += scaledFitnesses[fitIndex];
      }
      final float meanFitness = sumOfFitnesses / scaledFitnesses.length;
      System.out.println("Mean Scaled Fitness: " + meanFitness);
      float stddevOfFitnesses = 0.0f;
      for (int fitIndex = 0; fitIndex < scaledFitnesses.length; ++fitIndex) {
        stddevOfFitnesses += Math.pow((scaledFitnesses[fitIndex] - meanFitness), 2);
      }
      stddevOfFitnesses = (float) Math.sqrt(stddevOfFitnesses / scaledFitnesses.length);
      System.out.println("Std. Dev.: " + stddevOfFitnesses);
      
      // select the fittest individuals from the current generation
      ArrayList<String> selections = sim.selectFittest(scaledFitnesses);
      // Add the two elite individuals to the new population
      ArrayList<String> population = new ArrayList<String>(scaledFitnesses.length);
      population.add(selections.remove(0));
      population.add(selections.remove(0));
      // pair off and mate the fittest individuals to create the next generation
      System.out.println("Mating Selection Count: " + selections.size());
      final Random rand = new Random();
      while (selections.size() > 0) {
        String[] children = HprfRfHammer2DGeneticSolver.crossover(
          selections.remove((int) Math.round((selections.size() - 1) * rand.nextFloat())),
          selections.remove((int) Math.round((selections.size() - 1) * rand.nextFloat())));
        for (int childIndex = 0; childIndex < 2; ++childIndex) {
          population.add(HprfRfHammer2DGeneticSolver.mutate(children[childIndex], MUTATION_RATE));
        }
      }

      // add the next generation to the simulation
      sim.clearParams();
      Iterator<String> iter = population.iterator();
      while (iter.hasNext()) {
        sim.addParamsFromString(iter.next());
      }

      // run the simulation and generate a fitness for each individual
      fitnesses = sim.evaluate(trueData);
      scaledFitnesses = scaleFitnesses(fitnesses);

      System.out.println("### Generation " + generationIndex + " Fitnesses ###");
      final float lastAvgFitness = avgFitness;
      avgFitness = 0.0f;
      for (int fitIndex = 0; fitIndex < fitnesses.length; ++fitIndex) {
        avgFitness += fitnesses[fitIndex];

        if (scaledFitnesses[fitIndex] > scaledFitnesses[maxFitnessIndex]) {
          maxFitnessIndex = fitIndex;
        }
        System.out.println("Individual " + fitIndex + ": " + scaledFitnesses[fitIndex]);
        System.out.println(sim.getParamsAsString(fitIndex));
      }
      avgFitness /= fitnesses.length;
      System.out.println("Fittest Individual: " + maxFitnessIndex + " (" + scaledFitnesses[maxFitnessIndex] + ")");
      System.out.println(sim.getParamsAsString(maxFitnessIndex));

      final float deltaAvgFitness = Math.abs(avgFitness - lastAvgFitness);
      System.out.println("Change in Average Fitness: " + deltaAvgFitness);
      if (deltaAvgFitness < RELATIVE_TOLERANCE) {
        System.out.println("WARNING! Run aborted due to the average relative fitness being below " +
                           RELATIVE_TOLERANCE + ".");
      }
    }  // generation loop

    System.out.println("### Final Results ###");
    for (int individualIndex = 0; individualIndex < POPULATION_SIZE; ++individualIndex) {
      final StringBuffer line = new StringBuffer();
      line.append(individualIndex).append(" -- ").append(fitnesses[individualIndex])
          .append(" -- ");
      final float[] parameters = sim.getParams(individualIndex);
      for (int paramIndex = 0; paramIndex < parameters.length; ++paramIndex) {
        line.append(parameters[paramIndex]).append(" ");
      }
      System.out.println(line.toString());
    }

    float[][][] data = sim.getData();
    System.out.print("data = np.array((");
    for (int plotIndex = 0; plotIndex < 3; ++plotIndex) {
      if (plotIndex > 0) {
        System.out.println(",");
      }
      System.out.print("(");
      for (int sampleIndex = 0; sampleIndex < data[0][0].length; ++sampleIndex) {
        if (sampleIndex > 0) {
          System.out.print(",");
        }
        System.out.print(data[maxFitnessIndex][plotIndex][sampleIndex]);
      }
      System.out.print(")");
    }
    System.out.println("))");
  }

  public static void main(String[] args) {
    for (int argIndex = 0; argIndex < args.length; ++argIndex) {
      System.out.println("Arg " + argIndex + ": " + args[argIndex]);
    }
    if (args[0].equals("test")) {
      String[] testargs = new String[args.length-2];
      System.arraycopy(args, 2, testargs, 0, args.length-2);
      if (args[1].equals("main")) {
        System.out.println("Running testMain...");
        testMain(testargs);
      } else if (args[1].equals("fitness")) {
        System.out.println("Running testFitness...");
        testFitness(testargs);
      }
    } else {
      main0(args);
    }
  }
/**************** TODO ********************
x 0) Implement CSV data loading
x 1) Determine what settings will be varied.
    r_S0, r_S5, z_S1, rho_SS, rho_Cu, K_SS, K_Cu, G_SS, G_Cu
x 2) Create methods for [de]serializing variable settings to/from a string.
x 3) Create methods for mutating and mating the "DNA" string
Genetic Algorithm:
    x i)    Create a starter model and save it.
    ii)   HprfRfHammer2DGeneticSolver.load() N copies of the starter model
    iii)  Mutate the N models
    iv)   Build and run the N models, determining the fitness of each
    v)    Take the fittest P models. Stop if any of their fitness levels are above some threshold F. 
    vii)  Mate the P fittest models to generate N new models
    viii) Go to step (iii)
*/
}
