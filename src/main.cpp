// Qt
#include <QApplication>
#include <QDesktopWidget>
#include <QThread>

// sph
#include "sph.h"
#include "sphconfig.h"
#include "visualization.h"
#include "widget.h"


int main(int argc, char *argv[])
{
   QApplication a(argc, argv);

//   int width = a.desktop()->width();
//   int height = a.desktop()->height();

   SPH sph;

   // run simulation whitout visualization
   if (argc > 1 && std::string(argv[1]) == "r") 
   {
      sph.start();
      sph.wait();
      return 0;
   }

   Widget w;
   w.getVisualization()->setSph(&sph);
   w.show();
   // w.move(width/4, height/4);

   // init config widget
   w.Config()->setSph(&sph);
   w.Config()->readValuesFromSimulation();

   a.connect(
      &sph,
      SIGNAL(updateElapsed(int, int, int, int, int, int)),
      &w,
      SLOT(updateElapsedSph(int, int, int, int, int, int)),
      Qt::QueuedConnection
   );

   sph.start();
   // sph.pauseResume();

   a.connect(
      &w,
      SIGNAL(startClicked()),
      &sph,
      SLOT(pauseResume())
   );

   a.connect(
      &w,
      SIGNAL(shutDownClicked()),
      &sph,
      SLOT(stopSimulation())
   );

   a.connect(
      &sph,
      SIGNAL(finished()),
      &a,
      SLOT(quit())
   );

   return a.exec();
}
