import java.awt.Color;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
import javax.imageio.ImageIO;
import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ps.PSGraphics2D;

class CreateLogo
{
  public static void main(String[] paramArrayOfString)
  {
    CreateLogo localCreateLogo = new CreateLogo();

    String str1 = paramArrayOfString[0] + "project.txt";
    String str2 = paramArrayOfString[0] + "EM_project.txt";

    String str3 = paramArrayOfString[1].toUpperCase();
    int i;
    if (str3.equals("DNA"))
      i = 0;
    else {
      i = 1;
    }
    boolean bool = Boolean.parseBoolean(paramArrayOfString[4].toUpperCase());
    System.out.println("Unique?: " + bool);

    System.out.println("Loading data for: " + str1);
    System.out.println("Generating logos:");

    File localFile1 = new File(str1);

    List localList = ProfileReader.readAsProfiles(localFile1, 0.01D, bool, i);

    for (int j = 0; j < localList.size(); j++) {
      BindingProfile localObject = (BindingProfile)localList.get(j);
      localList.size();
      localCreateLogo.saveALogo((BindingProfile)localObject, j, paramArrayOfString[2], paramArrayOfString[3]);
    }

    System.out.println("Loading data for: " + str2);
    System.out.println("Generating logos:");
    
    File localFile2 = new File(str2);

    Object localObject = ProfileReader.readAsProfiles(localFile2, 0.01D, bool, i);

    for (int k = 0; k < ((List)localObject).size(); k++) {
      BindingProfile localBindingProfile = (BindingProfile)((List)localObject).get(k);
      ((List)localObject).size();
      localCreateLogo.saveALogo(localBindingProfile, k, paramArrayOfString[2], paramArrayOfString[3]);
    }

    System.out.println("Done.");
  }

  public void saveALogo(BindingProfile paramBindingProfile, int paramInt, String paramString1, String paramString2) {
    try {
      String str1 = paramBindingProfile.getName();
    
      System.out.println("\t" + str1);
      ProfileSequenceLogo localProfileSequenceLogo = new ProfileSequenceLogo(paramBindingProfile, 180);
	  
      String str2 = paramString2 + "/logos/" + str1;

      //int i = localProfileSequenceLogo.getFigureWidth();
      int i = localProfileSequenceLogo.getLogoWidth();
      int j = localProfileSequenceLogo.getLogoHeight();
      
      PSGraphics2D localPSGraphics2D = null;

      
      Object localObject;
      if (paramString1.equals("png")) {
        str2 = str2 + ".png";
        localObject = localProfileSequenceLogo.drawSequenceLogo();
        ImageIO.write((RenderedImage)localObject, "png", new File(str2));
      } else if (paramString1.equals("eps")) {
        str2 = str2 + ".eps";
        localPSGraphics2D = new PSGraphics2D(new File(str2), new Dimension(i, j));
        localPSGraphics2D.startExport();
        localPSGraphics2D.setBackground(Color.WHITE);
        localPSGraphics2D.clearRect(0, 0, i, j);
        localPSGraphics2D.setColor(Color.BLACK);
	
        localProfileSequenceLogo.drawSequenceLogo(localPSGraphics2D, new Point(i / 2, j / 2));
        localPSGraphics2D.endExport();
      } else if (paramString1.equals("pdf")) {

	  localObject = new FileOutputStream(str2 + ".pdf");
	  localProfileSequenceLogo.saveAsPDF((OutputStream)localObject, null);
	  
      } else if (paramString1.equals("pdfpng"))
      {
        localObject = new FileOutputStream(str2 + ".pdf");
        localProfileSequenceLogo.saveAsPDF((OutputStream)localObject, null);

        BufferedImage localBufferedImage = localProfileSequenceLogo.drawSequenceLogo();
        ImageIO.write(localBufferedImage, "png", new File(str2 + ".png"));
      } else {
        System.out.println("Unknown File Format: " + paramString1);
        System.exit(3);
      }
    }
    catch (Exception localException) {
      System.out.println("Exception: " + localException);
    }
  }
}
