#include "taskcontrolclient.h"
#include "ui_taskcontrolclient.h"

taskControlClient::taskControlClient(QWidget *parent) : QMainWindow(parent), ui(new Ui::taskControlClient) {
    ui->setupUi(this);
    my_socket = new QTcpSocket;

    maxlength=200;
    procNum = 0;
    setServerIPaddress(127,0,0,1,3333);
    QObject::connect(my_socket,SIGNAL(connected()),this,SLOT(on_connectedToTCPSocket()));
    QObject::connect(my_socket,SIGNAL(disconnected()),this,SLOT(on_disconnectedFromTCPSocket()));
    QObject::connect(my_socket,SIGNAL(error(QAbstractSocket::SocketError)),this,SLOT(on_TCPsocketConnectionError()));
    QObject::connect(my_socket,SIGNAL(readyRead()),this,SLOT(on_TCPsocketReadyToRead()));


    tab2grid =  new QGridLayout(ui->tab_2);
//    tab2grid->setSpacing(6);
    tab2grid->setObjectName(QString::fromUtf8("tab2grid"));
//    tab2grid->addLayout(tab2grid, 0, 0, 1, 2);
    tab2grid->setMargin(0);
    tab2grid->setHorizontalSpacing(0);
    tab2grid->setVerticalSpacing(0);
//    tab2grid->se

//    test = new progressBar(1,12,123,21,ui->tab_2);
//    test = new progressBar(1,12,123,21,ui->tab_2);
//    test->setGeometry(0,20,100,40);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    test->setGeometry(100,20,200,40);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    tab2grid->addWidget(test);


}

taskControlClient::~taskControlClient()
{
    delete ui;
    delete my_socket;
    for (int i=0;i<node.count();i++) delete run[i]; run.clear();
    delete tab2grid;
}



void taskControlClient::setServerIPaddress(int ip1, int ip2, int ip3, int ip4, int port) {
    ui->spinBox_RTserverIPaddress1->setValue(ip1);
    ui->spinBox_RTserverIPaddress2->setValue(ip2);
    ui->spinBox_RTserverIPaddress3->setValue(ip3);
    ui->spinBox_RTserverIPaddress4->setValue(ip4);
    ui->spinBox_RTserverPort->setValue(port);
}


void taskControlClient::server_changed() {
    if (ui->comboBox_serverSelection->currentText()=="cosmos") {       setServerIPaddress(158,75,6,55,3333); }
    if (ui->comboBox_serverSelection->currentText()=="szajcomp") {        setServerIPaddress(127,0,0,1,3333); }
}

void taskControlClient::on_pushButton_checkStatus_clicked() {
    printf("connecting to server\n");
    my_socket->connectToHost( getIPaddr(), (quint16)(ui->spinBox_RTserverPort->value()) );

    node.clear();
    stSim.clear();
    enSim.clear();
    curSim.clear();

    checkRunStatus();

}

void taskControlClient::on_connectedToTCPSocket() {
    printf("Connected to the server.\n");
    //        ui->pushButton_disconectFromServer->setEnabled(true);
    ui->pushButton_checkStatus->setEnabled(false);
    //    startDataSend();
    //        RTconnected=true;
}

void taskControlClient::on_disconnectedFromTCPSocket() {
    printf("Disconnected from the server.\n");
    ui->pushButton_checkStatus->setEnabled(true);
    //        RTconnected=false;
    printData();
    updateGUI();
}

void taskControlClient::on_TCPsocketConnectionError() {
    printf("Error connecting to the server. %s\n",getIPaddr().toStdString().c_str());
    ui->pushButton_checkStatus->setEnabled(true);
    //        RTconnected=false;

}

const QString taskControlClient::getIPaddr() const {
    QString s,tmps;
    s=tmps.setNum(ui->spinBox_RTserverIPaddress1->value());
    s+="."+tmps.setNum(ui->spinBox_RTserverIPaddress2->value())+".";
    s+=tmps.setNum(ui->spinBox_RTserverIPaddress3->value())+".";
    s+=tmps.setNum(ui->spinBox_RTserverIPaddress4->value());
//    printf("ip addr is: %s\n",s.toStdString().c_str());
    return s;
}

void taskControlClient::checkRunStatus() {
    QString s;
    char tmpch[maxlength];
    QString cmd="run status\r\n";
    qint64 data=my_socket->write(cmd.toStdString().c_str());

}

void taskControlClient::on_TCPsocketReadyToRead() {
    char tmpch[maxlength];
    qint64 data;
    QString s;
    QByteArray qba;
    printf("printing answer from the server\n");
    do {
//        data=my_socket->readLine(tmpch,maxlength);
        qba=my_socket->readLine();
//        if (data>0) {
        if (qba.count()>0) {
//            printf("*%s* %i,bytes left:%i\n",tmpch,int(data),int(my_socket->bytesAvailable()));
            printf("%s, read:%i bytes, bytes left:%i\n",qba.constData(),int(qba.size()),int(my_socket->bytesAvailable()));
//            s=tmpch;
            s=qba.constData();
            extractData(s);
        }
        else {
            printf("read no data from the socket\n");
        }
        //        data=my_socket->readLine(tmpch,maxlength);
//    }while (data>0);
    } while (my_socket->bytesAvailable()>0);

}

void taskControlClient::extractData(QString& data) {
    int id, st, en, cur, num;
    sscanf(data.toStdString().c_str(),"%i %i %i %i %i",&id, &num, &st, &en, &cur);

    node.append(id);
    stSim.append(st);
    enSim.append(en);
    curSim.append(cur);
    procNum=node.count();
}


void taskControlClient::printData() const {
    for (int i=0;i<node.count();i++) {
        printf("%i %i %i %i\n",node.value(i),stSim.value(i),enSim.value(i),curSim.value(i));
    }

}


void taskControlClient::updateGUI() {
//ui->tab_2->
//    printf("node count %i\n",node.count()); exit(0);
    for (int i=0;i<run.count();i++) delete run[i]; run.clear();
//exit(0);
    for (int i=0;i<procNum;i++) {
        printf("adding progress bar\n");
        run.append(new progressBar(i,stSim[i],enSim[i],curSim[i],ui->tab_2));
//        run.last()->setObjectName(QString::fromUtf8("progressBar"));
//        run.last()->setGeometry(QRect(130, 80, 113, 29));
        tab2grid->addWidget(run.last());
        run.last()->show();
//                //    test = new progressBar(1,12,123,21,ui->tab_2);
//    test = new progressBar(2,12,123,21,ui->tab_2);
//    tab2grid->addWidget(test);
//    test->setGeometry(0,20,300,40);
//    test->show();

    }
}
