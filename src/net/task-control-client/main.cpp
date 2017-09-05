#include <QtGui/QApplication>
#include "taskcontrolclient.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    taskControlClient w;
    w.show();
    return a.exec();
}
