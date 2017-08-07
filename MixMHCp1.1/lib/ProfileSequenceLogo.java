import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.DefaultFontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfTemplate;
import com.lowagie.text.pdf.PdfWriter;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.font.LineMetrics;
import java.awt.font.TextLayout;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.gui.SymbolStyle;
import org.biojava.bio.gui.TextLogoPainter;

public class ProfileSequenceLogo
{
    private BindingProfile profile = null;
    
    private int sequenceLogoStartIndex = 1;
    private boolean trimLogo = false;
    private double trimLogoPercentage = 0.1D;
    private int logoHeight = 240;
    private int figOffset=40;       //This is important to display the y axis (number of bits).

  private LogoSizeSpecification logoSpec = null;

  public ProfileSequenceLogo(BindingProfile profile, int logoHeight)
  {
    this.profile = profile;
    this.trimLogo = false;
    this.logoSpec = new LogoSizeSpecification(this, profile, logoHeight);
  }

  public ProfileSequenceLogo(BindingProfile profile, double trimLogoPercentage, int logoHeight)
  {
    this.profile = profile;
    this.trimLogoPercentage = trimLogoPercentage;
    this.trimLogo = true;
    this.logoHeight = logoHeight;
    this.logoSpec = new LogoSizeSpecification(this, profile, this.logoHeight);
  }

  public void sequenceLogoSetStartIndex(int startIndex)
  {
    this.sequenceLogoStartIndex = startIndex;
  }

  public boolean isTrimLogo()
  {
    return this.trimLogo;
  }

  public double getTrimLogoPercentage()
  {
    return this.trimLogoPercentage;
  }

  public LogoSizeSpecification getDetailedLogoSizeSpecification()
  {
    return this.logoSpec;
  }

  public int getLogoHeight()
  {
    return this.logoHeight;
  }

  public int getLogoWidth()
  {
    return (int)this.logoSpec.logoWidth;
  }
  public int getFigureWidth()
  {
    return (int)this.logoSpec.logoWidth + figOffset;
  }

  public BufferedImage drawSequenceLogo()
  {
    BufferedImage bi = new BufferedImage((int)this.logoSpec.logoWidth + figOffset, (int)this.logoSpec.logoHeight, 1);
    Graphics2D graphics = bi.createGraphics();
    graphics.setBackground(Color.WHITE);
    graphics.setColor(Color.WHITE);
    graphics.clearRect(0, 0, (int)this.logoSpec.logoWidth + figOffset, (int)this.logoSpec.logoHeight);
    drawSequenceLogo(graphics, null);
    return bi;
  }

  public BufferedImage drawSequenceLogo(SymbolStyle symbolColorStyle)
  {
    BufferedImage bi = new BufferedImage((int)this.logoSpec.logoWidth + figOffset, (int)this.logoSpec.logoHeight, 1);
    Graphics2D graphics = bi.createGraphics();
    graphics.setBackground(Color.WHITE);
    graphics.setColor(Color.WHITE);
    graphics.clearRect(0, 0, (int)this.logoSpec.logoWidth + figOffset, (int)this.logoSpec.logoHeight);
    drawSequenceLogo(graphics, null, symbolColorStyle);
    return bi;
  }

  public void saveAsPDF(OutputStream out, SymbolStyle symbolColorStyle) throws IOException
  {
    Rectangle pagesize = new Rectangle((float)this.logoSpec.logoWidth + figOffset, (float)this.logoSpec.logoHeight);
    Document document = new Document(pagesize, 50.0F, 50.0F, 50.0F, 50.0F);
    try
    {
      PdfWriter writer = PdfWriter.getInstance(document, out);
      document.open();
      PdfContentByte cb = writer.getDirectContent();
      PdfTemplate tp = cb.createTemplate((float)this.logoSpec.logoWidth + figOffset, (float)this.logoSpec.logoHeight);

      Graphics2D graphics = tp.createGraphics((float)this.logoSpec.logoWidth + figOffset, (float)this.logoSpec.logoHeight, new DefaultFontMapper());
      graphics.setBackground(Color.WHITE);
      graphics.setColor(Color.WHITE);
      graphics.clearRect(0, 0, (int)this.logoSpec.logoWidth + figOffset, (int)this.logoSpec.logoHeight);
      if (symbolColorStyle == null)
        drawSequenceLogo(graphics, null);
      else
        drawSequenceLogo(graphics, null, symbolColorStyle);
      graphics.dispose();

      cb.addTemplate(tp, 0.0F, 0.0F);
    }
    catch (DocumentException de) {
      System.err.println(de.getMessage());
    }
    document.close();
  }

  public void drawSequenceLogo(Graphics2D g, Point center)
  {
    if (this.profile.getType() == 1)
    {
      drawSequenceLogo(g, center, new WebLogoProteinStyle());
    }
    else
    {
      drawSequenceLogo(g, center, new WebLogoDNAStyle());
    }
  }

  public void drawSequenceLogo(Graphics2D gr, Point center, SymbolStyle symbolColorStyle)
  {
    int logoWidth = (int)this.logoSpec.logoWidth + figOffset;
    int logoHeight = (int)this.logoSpec.logoHeight;
    int columnWidth = (int)this.logoSpec.columnWidth;
    int columnHeight = (int)this.logoSpec.columnHeight;

    WeightMatrix wm = this.profile.getWeightMatrix();

    JPanel logoPanel = new JPanel(new GridLayout(1, this.logoSpec.numberOfColumns));
    Dimension logoPanelSize = new Dimension(columnWidth * this.logoSpec.numberOfColumns, columnHeight);
    logoPanel.setPreferredSize(logoPanelSize);
    logoPanel.setOpaque(false);
    RenderingHints hints = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
    try
    {
      for (int pos = this.logoSpec.minColumnIndex; pos <= this.logoSpec.maxColumnIndex; pos++) {
        Distribution dist = wm.getColumn(pos);
        DistributionLogo dl = new DistributionLogo();
        dl.setRenderingHints(hints);
        dl.setOpaque(false);

        dl.setDistribution(dist);
        dl.setPreferredSize(new Dimension(columnWidth, columnHeight));
        dl.setLogoPainter(new TextLogoPainter());
        dl.setStyle(symbolColorStyle);
        logoPanel.add(dl); //This draws each letter
      }
    } catch (BioException ex) {
      throw new BioError(ex);
    }

    ProfileSequenceLogo.LogoAxisPanel axisPanel = new ProfileSequenceLogo.LogoAxisPanel();
    axisPanel.setLogoSizeSpecification(this.logoSpec);
   
    Dimension axisPanelSize = new Dimension(logoWidth, logoHeight);
    axisPanel.setPreferredSize(axisPanelSize);
    axisPanel.setOpaque(false);
    
    JPanel panel = new JPanel();
    panel.setOpaque(false);
    panel.setPreferredSize(axisPanelSize);
    panel.setLayout(null);

    axisPanel.setBounds(0, 0, (int)axisPanelSize.getWidth(), (int)axisPanelSize.getHeight());
    panel.add(axisPanel); //This is for the x axis

    logoPanel.setBounds((int)(axisPanelSize.getWidth() - logoPanelSize.getWidth()), 0, (int)logoPanelSize.getWidth(), (int)logoPanelSize.getHeight());
    panel.add(logoPanel); //This is for the logo

    AffineTransform at = null;
    if (center != null) {
	 at = gr.getTransform();
	 gr.translate(center.getX() - logoWidth / 2, center.getY() - logoHeight / 2);
    }
    
    JFrame frame = new JFrame();
    frame.getContentPane().add(panel);
    frame.pack();
    panel.print(gr); //This is to print the logo+axes
    frame.dispose();

    if (at != null)
      gr.setTransform(at);
  }

  private void drawRotatedString(String string, Graphics2D g, int x, int y)
  {
    TextLayout layout = new TextLayout(string, g.getFont(), g.getFontRenderContext());

    AffineTransform startMatrix = g.getTransform();
    g.translate(x, y);
    g.rotate(-1.570796326794897D);
    layout.draw(g, -layout.getAdvance() / 2.0F, 0.0F);

    g.setTransform(startMatrix);
  }

  private class LogoAxisPanel extends JComponent
  {
      private LogoSizeSpecification logoSpec = null;
      private String title = null;
      

    private LogoAxisPanel() {  } 
    public void setLogoSizeSpecification(LogoSizeSpecification logoSpec) { this.logoSpec = logoSpec; }

    public void setTitle(String title)
    {
      this.title = title;
    }

    //Required for drawing the axes
    protected void paintComponent(Graphics gOrig) {
      Graphics2D g = (Graphics2D)gOrig;

      int logoWidth = (int)this.logoSpec.logoWidth;
      int logoHeight = (int)this.logoSpec.logoHeight;
      int padLeft = (int)this.logoSpec.padLeft + figOffset;  //This is the space to the left of the logo, which is useful to display the scale
      int padBottom = (int)this.logoSpec.padBottom;
      int columnWidth = (int)this.logoSpec.columnWidth;
      int columnHeight = (int)this.logoSpec.columnHeight;

      Dimension panelSize = new Dimension(logoWidth, logoHeight);

      g.setColor(Color.BLACK);
      g.setFont(new Font("Dialog", 1, (int)this.logoSpec.fontSize));
      int yAxisOffset = 3;
      int yAxisTickLength = (int)this.logoSpec.yAxisTickLength;
      int tickLabelPadding = 2;
      int labelWidth = 0;
      g.drawLine(padLeft - yAxisOffset, 0, padLeft - yAxisOffset, columnHeight);
      double maxBits = Math.log(this.logoSpec.sequenceAlphabetSize) / Math.log(2.0D);

      //This is drawing the marks on the y axis
      for (int i = 0; i < maxBits; i++) {
        int y = (int)(columnHeight - i * columnHeight / maxBits);
        g.drawLine(padLeft - yAxisTickLength, y, padLeft - yAxisOffset, y);
        FontMetrics f = g.getFontMetrics();
        String tickLabel = Integer.toString(i);
        int tickLabelXOffset = (int)f.getStringBounds(tickLabel, g).getWidth();
        labelWidth = tickLabelXOffset;
        LineMetrics l = f.getLineMetrics(tickLabel, g);
        int tickLabelYOffset = (int)Math.floor(l.getAscent() / 2.0F);
	g.drawString(tickLabel, padLeft - yAxisTickLength - tickLabelXOffset - tickLabelPadding, y + tickLabelYOffset);
      }

      int columnLabelPadding = 2;
      //String yAxisLabel = "bits";
      //int yAxisLabelXPos = padLeft - yAxisTickLength - labelWidth - tickLabelPadding * 2;
      //int yAxisLabelYPos = columnHeight / 2;
      //ProfileSequenceLogo.this.drawRotatedString(yAxisLabel, g, yAxisLabelXPos, yAxisLabelYPos);
      
      //Required for dawing the position labels on the x axis
      for (int i = this.logoSpec.minColumnIndex; i <= this.logoSpec.maxColumnIndex; i++) {
        String columnLabel = Integer.toString(ProfileSequenceLogo.this.sequenceLogoStartIndex + i);
        int x = padLeft + columnWidth / 2 + (i - this.logoSpec.minColumnIndex) * columnWidth;
        int y = columnHeight + columnLabelPadding;
        FontMetrics f = g.getFontMetrics();
        Rectangle2D stringBounds = f.getStringBounds(columnLabel, g);
        LineMetrics l = f.getLineMetrics(columnLabel, g);
        if (stringBounds.getWidth() >= columnWidth)
        {
          if (stringBounds.getWidth() >= padBottom)
          {
            AffineTransform startMatrix = g.getTransform();
            double scaleFactor = padBottom / stringBounds.getWidth();
            g.scale(scaleFactor, scaleFactor);

            y = (int)(y + stringBounds.getWidth() / 2.0D * scaleFactor);
            x = (int)(x + l.getAscent() / 2.0F * scaleFactor);
            ProfileSequenceLogo.this.drawRotatedString(columnLabel, g, (int)(x / scaleFactor), (int)(y / scaleFactor));

            g.setTransform(startMatrix);
          }
          else {
            y = (int)(y + stringBounds.getWidth() / 2.0D);
            x = (int)(x + l.getAscent() / 2.0F);
            ProfileSequenceLogo.this.drawRotatedString(columnLabel, g, x, y);
          }
        }
        else
        {
          y = (int)(y + l.getAscent());
          x = (int)(x - stringBounds.getWidth() / 2.0D);
          g.drawString(columnLabel, x, y);
        }

      }

      if (this.title != null) {
        int padTitle = 3;
        FontMetrics f = g.getFontMetrics();
        Rectangle2D stringBounds = f.getStringBounds(this.title, g);
        g.drawString(this.title, (int)(panelSize.getWidth() / 2.0D - stringBounds.getWidth() / 2.0D), (int)panelSize.getHeight() - padTitle);
      }
    }
  }
}
