#include <QtWidgets>
#include "MainWindow.h"
#include "SceneWidget.h"

MainWindow::MainWindow()
{
	// place a widget in the window
	QWidget *widget = new SceneWidget;
	setCentralWidget(widget);

	createActions();
	createMenus();

	QString message = tr("Example message.");
	statusBar()->showMessage(message);

	setWindowTitle(tr("Valley Winds"));
	setMinimumSize(800, 600);
	showMaximized();

	setAcceptDrops(true);
}

//#ifndef QT_NO_CONTEXTMENU
//	void MainWindow::contextMenuEvent(QContextMenuEvent *event)
//	{
//		QMenu menu(this);
//		menu.exec(event->globalPos());
//	}
//#endif // QT_NO_CONTEXTMENU


void MainWindow::about()
{
	//infoLabel->setText(tr("Invoked <b>Help|About</b>"));
	QMessageBox::about(this, tr("About Menu"),
		tr("The <b>Menu</b> example shows how to create "
			"menu-bar menus and context menus."));
}

void MainWindow::createActions()
{
	aboutAct = new QAction(tr("&About"), this);
	aboutAct->setStatusTip(tr("Show the application's About box"));
	connect(aboutAct, &QAction::triggered, this, &MainWindow::about);
}

void MainWindow::createMenus()
{
	helpMenu = menuBar()->addMenu(tr("&Help"));
	helpMenu->addAction(aboutAct);
}
