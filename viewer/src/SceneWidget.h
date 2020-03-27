#pragma once

#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>
#include <vtkGenericOpenGLRenderWindow.h>

class SceneWidget : public QVTKOpenGLNativeWidget
{
public:
	SceneWidget();

private:

	void CreateTestScene();

	vtkSmartPointer<vtkGenericOpenGLRenderWindow> mRenderWindow;
};