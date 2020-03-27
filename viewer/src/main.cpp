#include <QApplication.h>
#include "MainWindow.h"

// ---------------------------------------
int main(int argc, char* argv[])
{
	// create the application
	QApplication app(argc, argv);

	// create the main window
	MainWindow window;

	// show the window and run the application
	window.show();
	return app.exec();
}
