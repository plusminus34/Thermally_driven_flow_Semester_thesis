#pragma once

#include <QMainWindow>

class QAction;
class QMenu;
class QSplitter;
class SceneWidget;

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow();

	SceneWidget* GetSceneWidget() { return mSceneWidget; }



private slots:
	void about();

private:
	void createActions();
	void createMenus();

	// menus
	QMenu *helpMenu;
	QAction *aboutAct;

	// widgets
	SceneWidget* mSceneWidget;
};