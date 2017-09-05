#ifndef TASKCONTROLCLIENT_H
#define TASKCONTROLCLIENT_H

#include <QtGui/QMainWindow>
#include <QtGui/QLineEdit>
#include <QtNetwork/QTcpSocket>
#include <QtCore/QTimer>
//#include <QLineEdit>
#include <QGridLayout>
#include "progressbar.h"

namespace Ui
{
    class taskControlClient;
}

class taskControlClient : public QMainWindow
{
    Q_OBJECT

public:
    taskControlClient(QWidget *parent = 0);
    ~taskControlClient();

private:
    Ui::taskControlClient *ui;

    void checkRunStatus();
    void setServerIPaddress(int ip1, int ip2, int ip3, int ip4, int port);
    const QString getIPaddr() const;
    void extractData(QString& data);
    void printData() const;
    void updateGUI();

    qint64 maxlength;

    int procNum;
    QList<int> node, stSim, enSim, curSim;

    QTcpSocket *my_socket;
    QList<progressBar*> run;
    progressBar* test;

    QGridLayout *tab2grid;

private slots:

    void server_changed();
    void on_connectedToTCPSocket();
    void on_TCPsocketConnectionError();
    void on_disconnectedFromTCPSocket();
    void on_pushButton_checkStatus_clicked();
    void on_TCPsocketReadyToRead();


};

#endif // TASKCONTROLCLIENT_H
